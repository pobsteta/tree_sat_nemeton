# ==============================================================================
# TreeSatAI-Time-Series — Téléchargement du dataset depuis HuggingFace
#
# Dataset : IGNF/TreeSatAI-Time-Series
# URL     : https://huggingface.co/datasets/IGNF/TreeSatAI-Time-Series
# Licence : CC-BY-SA 4.0
#
# 50 381 patchs de 60×60 m avec séries temporelles S1+S2 sur 1 an
# 20 espèces européennes (15 genres) — Basse-Saxe, Allemagne
#
# Structure HDF5 (sentinel-ts) :
#   sen-2-data   : (T, 10, 6, 6) — S2 L2A réflectance 10 bandes
#   sen-1-asc-data : (T, 2, 6, 6) — S1 VV/VH orbite ascendante
#   sen-1-des-data : (T, 2, 6, 6) — S1 VV/VH orbite descendante
#   sen-2-masks  : (T, 2, 6, 6) — probabilité neige/nuage
#
# Labels : TreeSatBA_v9_60m_multi_labels.json
#   {"Fagus_sylvatica_6_154201_BI_NLF.tif": [["Fagus", 0.8], ["Quercus", 0.2]]}
# ==============================================================================

# HuggingFace API base URL pour le repo
HF_REPO_ID <- "IGNF/TreeSatAI-Time-Series"
HF_API_BASE <- "https://huggingface.co/api/datasets"
HF_RESOLVE_BASE <- "https://huggingface.co/datasets"

# Téléchargement robuste avec reprise pour gros fichiers
# Utilise curl -C - pour reprendre un téléchargement interrompu
.download_with_resume <- function(url, dest_file, file_name, size_mb,
                                   max_retries = 5) {
  use_curl <- nzchar(Sys.which("curl"))
  partial  <- paste0(dest_file, ".partial")

  for (attempt in seq_len(max_retries)) {
    if (attempt > 1) {
      wait <- min(2^attempt, 30)
      cli::cli_alert_warning("  Tentative {attempt}/{max_retries} dans {wait}s...")
      Sys.sleep(wait)
    }

    tryCatch({
      if (use_curl) {
        # curl avec -C - = reprise automatique, -L = suivre redirections
        # --retry = retries internes curl, --connect-timeout = timeout connexion
        res <- system2("curl", c(
          "-L", "-C", "-",
          "--retry", "3",
          "--retry-delay", "5",
          "--connect-timeout", "30",
          "--max-time", "0",
          "-o", shQuote(partial),
          shQuote(url)
        ), stdout = TRUE, stderr = TRUE)

        exit_code <- attr(res, "status")
        if (is.null(exit_code)) exit_code <- 0L

        if (exit_code == 0 && file.exists(partial)) {
          # Vérifier que le fichier partiel a une taille raisonnable
          actual_mb <- round(file.info(partial)$size / 1024 / 1024, 1)
          if (size_mb > 0 && actual_mb < size_mb * 0.95) {
            cli::cli_alert_warning("  Incomplet : {actual_mb}/{size_mb} Mo")
            next
          }
          file.rename(partial, dest_file)
          cli::cli_alert_success("  {file_name} OK ({actual_mb} Mo)")
          return(TRUE)
        } else {
          msg <- paste(tail(res, 3), collapse = " ")
          cli::cli_alert_warning("  curl exit {exit_code}: {msg}")
        }
      } else {
        # Fallback : download.file avec timeout augmenté
        old_timeout <- getOption("timeout")
        options(timeout = max(3600, size_mb * 2))
        on.exit(options(timeout = old_timeout), add = TRUE)
        download.file(url, dest_file, mode = "wb", quiet = FALSE)
        cli::cli_alert_success("  {file_name} OK")
        return(TRUE)
      }
    }, error = function(e) {
      cli::cli_alert_warning("  Erreur : {e$message}")
    })
  }

  # Toutes les tentatives échouées
  cli::cli_alert_danger("  Echec apr\u00e8s {max_retries} tentatives pour {file_name}")
  if (file.exists(partial)) {
    actual_mb <- round(file.info(partial)$size / 1024 / 1024, 1)
    cli::cli_alert_info("  Fichier partiel conserv\u00e9 ({actual_mb} Mo) : {partial}")
    cli::cli_alert_info("  Relancez download_treesatai_hf() pour reprendre")
  }
  return(FALSE)
}

# Extraction de zip avec progression — compatible Windows (ZIP64, chemins longs)
# R unzip() ne supporte pas les ZIP > 4 Go sur Windows.
# Stratégie : 7z > tar (bsdtar Windows 10+) > PowerShell > R unzip (fallback)
.unzip_with_progress <- function(zip_path, exdir) {
  zip_name <- basename(zip_path)
  zip_mb   <- round(file.info(zip_path)$size / 1024 / 1024, 1)
  is_large <- file.info(zip_path)$size > 500 * 1024 * 1024  # > 500 Mo

  dir.create(exdir, showWarnings = FALSE, recursive = TRUE)

  # Compter les fichiers dans le zip (via le listing R, ok même pour gros zips)
  n_files <- tryCatch({
    fl <- unzip(zip_path, list = TRUE)
    total_mb <- round(sum(fl$Length) / 1024 / 1024, 1)
    cli::cli_alert_info("Extraction : {zip_name} ({zip_mb} Mo) \u2192 {nrow(fl)} fichiers ({total_mb} Mo)")
    nrow(fl)
  }, error = function(e) {
    cli::cli_alert_info("Extraction : {zip_name} ({zip_mb} Mo)")
    0L
  })
  flush.console()

  # Pour les petits zips : extraction directe
  if (!is_large) {
    tryCatch({
      unzip(zip_path, exdir = exdir, overwrite = TRUE)
      n_out <- length(list.files(exdir, recursive = TRUE))
      cli::cli_alert_success("  {zip_name} extrait ({n_out} fichiers)")
      return(invisible(TRUE))
    }, error = function(e) {
      cli::cli_alert_danger("  Echec unzip : {e$message}")
      return(invisible(FALSE))
    })
  }

  # --- Gros zips : outil externe avec monitoring ---
  # Chemins absolus Windows-safe (guillemets, pas de short names)
  zip_win <- normalizePath(zip_path, mustWork = TRUE)
  exdir_win <- normalizePath(exdir, mustWork = TRUE)
  existing_before <- length(list.files(exdir, recursive = TRUE))

  # Choisir le meilleur outil disponible
  tool <- .find_unzip_tool()
  cli::cli_alert_info("  Outil d'extraction : {tool$name}")
  flush.console()

  # Construire la commande
  if (tool$name == "7z") {
    cmd <- tool$path
    args <- c("x", "-y", paste0("-o", exdir_win), zip_win)
  } else if (tool$name == "tar") {
    cmd <- tool$path
    args <- c("-xf", zip_win, "-C", exdir_win)
  } else if (tool$name == "powershell") {
    cmd <- tool$path
    args <- c("-NoProfile", "-Command",
      sprintf('Expand-Archive -Path "%s" -DestinationPath "%s" -Force',
              zip_win, exdir_win))
  } else {
    # Fallback R unzip (peut échouer sur ZIP64)
    cli::cli_alert_warning("  Aucun outil externe trouv\u00e9, utilisation de R unzip()...")
    cli::cli_alert_warning("  Pour les ZIP > 4 Go, installez 7-Zip : https://7-zip.org/")
    flush.console()
    tryCatch({
      unzip(zip_path, exdir = exdir, overwrite = TRUE)
      n_out <- length(list.files(exdir, recursive = TRUE))
      cli::cli_alert_success("  {zip_name} extrait ({n_out} fichiers)")
      return(invisible(TRUE))
    }, error = function(e) {
      cli::cli_alert_danger("  Echec : {e$message}")
      return(invisible(FALSE))
    })
  }

  # Lancer en arrière-plan et surveiller
  t_start <- Sys.time()
  log_file <- tempfile(fileext = ".log")

  system2(cmd, args, wait = FALSE, stdout = log_file, stderr = log_file)

  cli::cli_alert_info("  Extraction en cours...")
  flush.console()

  last_pct <- -1L
  stall_count <- 0L
  last_count <- existing_before

  repeat {
    Sys.sleep(5)

    current_count <- tryCatch(
      length(list.files(exdir, recursive = TRUE)),
      error = function(e) last_count
    )
    new_files <- current_count - existing_before
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    # Détecter si l'extraction a fini (plus de nouveaux fichiers + processus terminé)
    if (current_count == last_count) {
      stall_count <- stall_count + 1L
    } else {
      stall_count <- 0L
      last_count <- current_count
    }

    # Afficher la progression
    if (n_files > 0) {
      pct <- min(100L, round(100 * new_files / n_files))
    } else {
      pct <- 0L
    }

    if ((pct > last_pct || stall_count == 0) && new_files > 0 && elapsed > 0) {
      speed <- round(new_files / elapsed)
      if (speed > 0 && n_files > 0) {
        remaining <- max(0, round((n_files - new_files) / speed))
        mins <- remaining %/% 60
        secs <- remaining %% 60
        cli::cli_alert_info("  [{pct}%] {new_files}/{n_files} fichiers ({speed}/s, reste ~{mins}m{secs}s)")
      } else {
        cli::cli_alert_info("  {new_files} fichiers extraits...")
      }
      flush.console()
      last_pct <- pct
    }

    # Fin : plus de mouvement depuis 30s (6 checks) = processus terminé
    if (stall_count >= 6 && new_files > 0) break

    # Timeout de sécurité : 2h
    if (elapsed > 7200) {
      cli::cli_alert_danger("  Timeout apr\u00e8s 2h d'extraction")
      break
    }
  }

  # Vérifier les erreurs dans le log
  if (file.exists(log_file)) {
    log_content <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) "")
    errors <- grep("error|erreur|fail|cannot|impossible", log_content,
                    ignore.case = TRUE, value = TRUE)
    if (length(errors) > 0) {
      cli::cli_alert_warning("  Avertissements extraction :")
      for (err in utils::head(errors, 5)) cli::cli_alert_warning("    {err}")
    }
    unlink(log_file)
  }

  final_count <- length(list.files(exdir, recursive = TRUE)) - existing_before
  elapsed_total <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")))
  mins <- elapsed_total %/% 60
  secs <- elapsed_total %% 60
  cli::cli_alert_success("  {zip_name} extrait ({final_count} fichiers en {mins}m{secs}s)")
  invisible(TRUE)
}

# Trouver le meilleur outil d'extraction disponible
.find_unzip_tool <- function() {
  # 1. 7-Zip (le plus fiable pour gros ZIP sur Windows)
  if (.Platform$OS.type == "windows") {
    sz_paths <- c(
      Sys.which("7z"),
      "C:/Program Files/7-Zip/7z.exe",
      "C:/Program Files (x86)/7-Zip/7z.exe"
    )
    for (p in sz_paths) {
      if (nzchar(p) && file.exists(p)) return(list(name = "7z", path = p))
    }
  } else {
    sz <- Sys.which("7z")
    if (nzchar(sz)) return(list(name = "7z", path = sz))
  }

  # 2. tar (bsdtar sur Windows 10+, gère le zip)
  tar_path <- Sys.which("tar")
  if (nzchar(tar_path)) return(list(name = "tar", path = tar_path))

  # 3. PowerShell (Windows, lent mais fiable)
  if (.Platform$OS.type == "windows") {
    ps <- Sys.which("powershell")
    if (nzchar(ps)) return(list(name = "powershell", path = ps))
  }

  # 4. Aucun outil externe
  list(name = "r_unzip", path = NULL)
}

#' Télécharger le dataset TreeSatAI-Time-Series depuis HuggingFace
#'
#' Télécharge les données du dataset IGNF/TreeSatAI-Time-Series hébergé
#' sur HuggingFace. Par défaut, télécharge les labels, le split et les
#' séries temporelles sentinel-ts. Le dataset complet fait ~40 Go.
#'
#' @param dest_dir Répertoire de destination
#' @param components Composants à télécharger : "labels", "split", "sentinel-ts",
#'   "sentinel", "aerial", "geojson", ou "all"
#' @param overwrite Écraser les fichiers existants
#' @return Chemin vers le répertoire de données (invisible)
#' @export
#'
#' @examples
#' \dontrun{
#' # Télécharger labels + split (léger, ~10 Mo)
#' download_treesatai_hf(components = c("labels", "split"))
#'
#' # Télécharger tout le dataset sentinel time-series (~30 Go)
#' download_treesatai_hf(components = c("labels", "split", "sentinel-ts"))
#'
#' # Tout télécharger (~40 Go)
#' download_treesatai_hf(components = "all")
#' }
download_treesatai_hf <- function(dest_dir   = file.path(DATA_DIR, "treesatai"),
                                   components = c("labels", "split", "geojson"),
                                   overwrite  = FALSE) {

  init_project_dirs()
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  cli::cli_h1("TreeSatAI-Time-Series (IGNF) \u2014 T\u00e9l\u00e9chargement HuggingFace")
  cli::cli_text("")
  cli::cli_alert_info("Source : {.url https://huggingface.co/datasets/{HF_REPO_ID}}")
  cli::cli_alert_info("Licence : CC-BY-SA 4.0")
  cli::cli_alert_info("50 381 patchs \u00d7 20 esp\u00e8ces \u2014 S\u00e9ries temporelles S1+S2")
  cli::cli_text("")

  if ("all" %in% components) {
    components <- c("labels", "split", "geojson", "sentinel-ts", "sentinel", "aerial")
  }

  cli::cli_alert_info("Composants demand\u00e9s : {paste(components, collapse = ', ')}")
  cli::cli_text("")

  # --- Lister les fichiers du repo HuggingFace ---
  cli::cli_alert_info("Interrogation de l'API HuggingFace...")

  api_url <- paste0(HF_API_BASE, "/", HF_REPO_ID, "/tree/main")
  resp <- tryCatch(
    httr2::request(api_url) |>
      httr2::req_headers("Accept" = "application/json") |>
      httr2::req_perform(),
    error = function(e) {
      cli::cli_alert_danger("Impossible de contacter HuggingFace : {e$message}")
      return(NULL)
    }
  )

  if (is.null(resp)) return(invisible(dest_dir))

  repo_files <- httr2::resp_body_json(resp)
  cli::cli_alert_success("{length(repo_files)} fichiers/dossiers trouv\u00e9s dans le repo")

  # --- Téléchargement de chaque composant ---
  for (component in components) {
    cli::cli_h2("T\u00e9l\u00e9chargement : {component}")

    # Trouver les fichiers correspondants
    matching <- Filter(function(f) {
      grepl(component, f$path, ignore.case = TRUE)
    }, repo_files)

    if (length(matching) == 0) {
      # Essayer de lister le sous-dossier
      subdir_url <- paste0(HF_API_BASE, "/", HF_REPO_ID, "/tree/main/", component)
      sub_resp <- tryCatch(
        httr2::request(subdir_url) |>
          httr2::req_headers("Accept" = "application/json") |>
          httr2::req_perform(),
        error = function(e) NULL
      )

      if (!is.null(sub_resp)) {
        matching <- httr2::resp_body_json(sub_resp)
      }
    }

    if (length(matching) == 0) {
      cli::cli_alert_warning("Aucun fichier trouv\u00e9 pour '{component}'")
      next
    }

    comp_dir <- file.path(dest_dir, component)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    for (file_info in matching) {
      file_path <- if (!is.null(file_info$path)) file_info$path else file_info$rfilename
      if (is.null(file_path)) next

      file_size <- file_info$size %||% 0
      file_name <- basename(file_path)
      dest_file <- file.path(dest_dir, file_path)

      dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)

      if (file.exists(dest_file) && !overwrite) {
        # Vérifier que le fichier n'est pas tronqué
        actual_size <- file.info(dest_file)$size
        if (file_size > 0 && actual_size >= file_size * 0.95) {
          cli::cli_alert_info("  D\u00e9j\u00e0 pr\u00e9sent : {file_name}")
          next
        } else if (file_size > 0) {
          cli::cli_alert_warning("  {file_name} incomplet ({round(actual_size/1024/1024,1)}/{round(file_size/1024/1024,1)} Mo), reprise...")
        } else {
          cli::cli_alert_info("  D\u00e9j\u00e0 pr\u00e9sent : {file_name}")
          next
        }
      }

      # URL de téléchargement direct
      dl_url <- paste0(HF_RESOLVE_BASE, "/", HF_REPO_ID, "/resolve/main/", file_path)

      size_mb <- round(file_size / 1024 / 1024, 1)
      cli::cli_alert_info("  T\u00e9l\u00e9chargement : {file_name} ({size_mb} Mo)")

      # Pour les gros fichiers (> 100 Mo), utiliser curl avec reprise
      dl_ok <- .download_with_resume(dl_url, dest_file, file_name, size_mb)
    }
  }

  # --- Post-traitement : extraction des zips si présents ---
  zip_files <- list.files(dest_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
  if (length(zip_files) > 0) {
    cli::cli_h2("Extraction des archives")
    for (zf in zip_files) {
      .unzip_with_progress(zf, dirname(zf))
    }
  }

  # --- Résumé ---
  cli::cli_h1("T\u00e9l\u00e9chargement termin\u00e9")
  cli::cli_text("")
  cli::cli_alert_success("Donn\u00e9es dans : {.path {dest_dir}}")

  all_files <- list.files(dest_dir, recursive = TRUE)
  cli::cli_alert_info("{length(all_files)} fichiers t\u00e9l\u00e9charg\u00e9s")

  # Vérifier les labels
  labels_file <- list.files(dest_dir, pattern = "multi_labels.*\\.json$",
                             recursive = TRUE, full.names = TRUE)
  if (length(labels_file) > 0) {
    labels <- jsonlite::fromJSON(labels_file[1])
    cli::cli_alert_success("{length(labels)} patchs labellis\u00e9s trouv\u00e9s")
    cli::cli_text("")
    cli::cli_text("Pour entra\u00eener un mod\u00e8le sur ces donn\u00e9es :")
    cli::cli_text("  {.code train_treesatai(data_path = \"{dest_dir}\")}")
  }

  invisible(dest_dir)
}

#' Charger les labels TreeSatAI depuis le JSON multi-labels
#'
#' @param labels_path Chemin vers TreeSatBA_v9_60m_multi_labels.json
#' @return data.frame avec colonnes : patch_id, species, genus, proportion
#' @export
load_treesatai_labels <- function(labels_path = NULL) {
  if (is.null(labels_path)) {
    labels_path <- list.files(
      file.path(DATA_DIR, "treesatai"),
      pattern = "multi_labels.*\\.json$",
      recursive = TRUE, full.names = TRUE
    )
    if (length(labels_path) == 0) {
      stop("Labels non trouv\u00e9s. Lancez d'abord download_treesatai_hf()")
    }
    labels_path <- labels_path[1]
  }

  cli::cli_alert_info("Chargement des labels : {labels_path}")
  raw <- jsonlite::fromJSON(labels_path, simplifyVector = FALSE)

  # Convertir en data.frame
  rows <- lapply(names(raw), function(patch_name) {
    # Extraire espèce dominante du nom de fichier
    # Format : Genus_species_ageclass_ID_dataset_source.tif
    parts <- strsplit(patch_name, "_")[[1]]
    genus <- parts[1]
    species_full <- paste(parts[1], parts[2], sep = "_")

    # Labels multi-genres avec proportions
    label_list <- raw[[patch_name]]
    do.call(rbind, lapply(label_list, function(lbl) {
      data.frame(
        patch_id   = patch_name,
        species    = species_full,
        genus      = lbl[[1]],
        proportion = lbl[[2]],
        stringsAsFactors = FALSE
      )
    }))
  })

  result <- do.call(rbind, rows)
  cli::cli_alert_success("{nrow(result)} \u00e9tiquettes charg\u00e9es ({length(unique(result$patch_id))} patchs)")
  result
}

#' Charger le split train/test TreeSatAI
#'
#' @param split_dir Répertoire contenant les fichiers .lst
#' @return Liste avec train et test (vecteurs de noms de patchs)
#' @export
load_treesatai_split <- function(split_dir = NULL) {
  if (is.null(split_dir)) {
    split_dir <- file.path(DATA_DIR, "treesatai", "split")
  }

  train_file <- list.files(split_dir, pattern = "train.*\\.lst$", full.names = TRUE)
  test_file  <- list.files(split_dir, pattern = "test.*\\.lst$", full.names = TRUE)

  train_patches <- if (length(train_file) > 0) readLines(train_file[1]) else character(0)
  test_patches  <- if (length(test_file) > 0)  readLines(test_file[1])  else character(0)

  cli::cli_alert_success("Split charg\u00e9 : {length(train_patches)} train / {length(test_patches)} test")

  list(train = train_patches, test = test_patches)
}

# %||% est défini dans 08_download_satellite.R
