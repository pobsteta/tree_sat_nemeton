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
    # Filtrer les zips déjà extraits (marqueur .extracted)
    needs_extract <- vapply(zip_files, function(zf) {
      marker <- paste0(zf, ".extracted")
      !file.exists(marker)
    }, logical(1))

    if (any(needs_extract)) {
      cli::cli_h2("Extraction des archives")
      for (zf in zip_files[needs_extract]) {
        ok <- .unzip_with_progress(zf, dirname(zf))
        if (isTRUE(ok)) {
          # Marquer comme extrait pour éviter de re-extraire
          writeLines(format(Sys.time()), paste0(zf, ".extracted"))
        }
      }
    }

    skipped <- zip_files[!needs_extract]
    if (length(skipped) > 0) {
      for (zf in skipped) {
        cli::cli_alert_info("D\u00e9j\u00e0 extrait : {basename(zf)}")
      }
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

# ==============================================================================
# Correspondance Genre → espèce (pour la classification au niveau genre)
# ==============================================================================

# Table de correspondance genre → nom français + code SPECIES
# Les genres sans correspondance SPECIES utilisent le code du genre le plus
# similaire pour la génération de profils synthétiques
GENUS_MAP <- list(
  Fagus        = list(french = "H\u00eatre",                code = 5),
  Quercus      = list(french = "Ch\u00eane p\u00e9doncul\u00e9", code = 1),
  Picea        = list(french = "\u00c9pic\u00e9a commun",   code = 13),
  Pinus        = list(french = "Pin sylvestre",             code = 16),
  Abies        = list(french = "Sapin pectin\u00e9",        code = 14),
  Larix        = list(french = "M\u00e9l\u00e8ze d'Europe", code = 20),
  Betula       = list(french = "Bouleau verruqueux",        code = 8),
  Acer         = list(french = "\u00c9rable sycomore",      code = 10),
  Fraxinus     = list(french = "Fr\u00eane commun",         code = 9),
  Carpinus     = list(french = "Charme",                    code = 7),
  Pseudotsuga  = list(french = "Douglas",                   code = 15),
  Populus      = list(french = "Peupliers",                 code = 11),
  Robinia      = list(french = "Robinier faux-acacia",      code = 12),
  Castanea     = list(french = "Ch\u00e2taignier",          code = 6),
  Alnus        = list(french = "Aulne",                     code = 8),   # similaire Bouleau
  Tilia        = list(french = "Tilleul",                   code = 10),  # similaire Érable
  Prunus       = list(french = "Prunus",                    code = 6),   # similaire Châtaignier
  Sorbus       = list(french = "Sorbier",                   code = 9),   # similaire Frêne
  Salix        = list(french = "Saule",                     code = 11),  # similaire Peuplier
  Ulmus        = list(french = "Orme",                      code = 7),   # similaire Charme
  # Classe non-arbre du dataset TreeSatAI original (coupe rase, vide forestier)
  Cleared      = list(french = "Coupe/Vide",                code = 21)
)

# Classes non forestières à exclure de l'entraînement
# Note : "Cleared" est CONSERVÉ car il fait partie des 15 genres TreeSatAI.
# C'est une classe importante pour détecter les coupes rases et les vides.
# Seules les classes vraiment hors-sujet sont exclues.
NON_TREE_CLASSES <- c("NonForest", "non_forest", "Water", "Urban")

# ==============================================================================
# Chargement complet des données TreeSatAI pour le pipeline
# ==============================================================================

#' Charger les données TreeSatAI et construire ts_long pour le pipeline
#'
#' Détecte le format TreeSatAI (labels JSON + patches HDF5 dans sentinel-ts),
#' lit les séries temporelles et construit un data.frame au format ts_long
#' compatible avec build_feature_matrix().
#'
#' @param data_path Chemin vers le répertoire TreeSatAI
#' @return Liste avec ts_long (data.frame), split (list train/test),
#'   genus_names (vecteur des genres uniques)
#' @export
load_treesatai_data <- function(data_path) {
  cli::cli_h2("Chargement des donn\u00e9es TreeSatAI")

  # --- 1. Labels ---
  labels_file <- list.files(data_path, pattern = "multi_labels.*\\.json$",
                             recursive = TRUE, full.names = TRUE)
  if (length(labels_file) == 0) {
    stop("Labels TreeSatAI non trouv\u00e9s dans ", data_path,
         ". Lancez d'abord download_treesatai_hf()")
  }
  labels <- load_treesatai_labels(labels_file[1])

  # Genre dominant par patch (proportion la plus élevée)
  dominant <- do.call(rbind, lapply(split(labels, labels$patch_id), function(df) {
    idx <- which.max(df$proportion)
    data.frame(
      patch_id     = df$patch_id[1],
      genus        = df$genus[idx],
      proportion   = df$proportion[idx],
      stringsAsFactors = FALSE
    )
  }))
  rownames(dominant) <- NULL

  # Mapper genre → nom français
  dominant$species_name <- vapply(dominant$genus, function(g) {
    info <- GENUS_MAP[[g]]
    if (!is.null(info)) info$french else g
  }, character(1))

  dominant$species_code <- vapply(dominant$genus, function(g) {
    info <- GENUS_MAP[[g]]
    if (!is.null(info) && !is.na(info$code)) as.integer(info$code) else 0L
  }, integer(1))

  # Filtrer les classes non forestières (Cleared, etc.)
  non_tree <- dominant$species_name %in% NON_TREE_CLASSES | dominant$genus %in% NON_TREE_CLASSES
  if (any(non_tree)) {
    n_removed <- sum(non_tree)
    removed_names <- paste(unique(dominant$species_name[non_tree]), collapse = ", ")
    cli::cli_alert_warning("Exclusion de {n_removed} patchs non forestiers ({removed_names})")
    dominant <- dominant[!non_tree, ]
  }

  genus_names <- sort(unique(dominant$species_name))
  cli::cli_alert_success("{nrow(dominant)} patchs, {length(genus_names)} genres : {paste(genus_names, collapse = ', ')}")

  # --- 2. Split train/test ---
  split_dir <- file.path(data_path, "split")
  split_info <- if (dir.exists(split_dir)) {
    load_treesatai_split(split_dir)
  } else {
    cli::cli_alert_warning("Split non trouv\u00e9, s\u00e9paration al\u00e9atoire 70/30")
    NULL
  }

  # --- 3. Recherche des données sentinel-ts (HDF5) ---
  ts_dir <- file.path(data_path, "sentinel-ts")
  ts_long <- NULL

  if (dir.exists(ts_dir)) {
    # Chercher les fichiers HDF5
    h5_files <- list.files(ts_dir, pattern = "\\.(h5|hdf5|hdf)$",
                            recursive = TRUE, full.names = TRUE)

    if (length(h5_files) > 0 && requireNamespace("hdf5r", quietly = TRUE)) {
      cli::cli_alert_info("Lecture de {length(h5_files)} patchs HDF5...")
      ts_long <- .read_treesatai_hdf5(h5_files, dominant)
    } else if (length(h5_files) > 0) {
      cli::cli_alert_warning(
        "Patchs HDF5 trouv\u00e9s mais {.pkg hdf5r} non install\u00e9. ",
        "install.packages('hdf5r') pour lire les donn\u00e9es r\u00e9elles."
      )
    }
  }

  # --- 4. Fallback : données synthétiques basées sur les labels réels ---
  if (is.null(ts_long)) {
    cli::cli_alert_info("G\u00e9n\u00e9ration de s\u00e9ries temporelles synth\u00e9tiques bas\u00e9es sur les labels r\u00e9els")
    ts_long <- .generate_ts_from_labels(dominant)
  }

  list(
    ts_long      = ts_long,
    split        = split_info,
    genus_names  = genus_names,
    dominant     = dominant
  )
}

# --- Lecture HDF5 des patches TreeSatAI ---
.read_treesatai_hdf5 <- function(h5_files, dominant) {
  # Dates cibles (73 pas de 5 jours sur 1 an)
  target_dates <- seq.Date(
    as.Date(TS_PARAMS$start_date),
    as.Date(TS_PARAMS$end_date),
    by = TS_PARAMS$target_interval_days
  )

  # Index des patchs par nom de fichier
  patch_names <- dominant$patch_id
  # Le nom HDF5 peut être Genus_species_age_ID_dataset_source.h5
  # Le label utilise .tif, le fichier est .h5
  h5_basenames <- tools::file_path_sans_ext(basename(h5_files))

  cli::cli_progress_bar("Lecture HDF5", total = length(h5_files))

  all_rows <- list()

  for (i in seq_along(h5_files)) {
    cli::cli_progress_update()

    h5_name <- h5_basenames[i]

    # Chercher le patch dans les labels (nom .tif ou sans extension)
    match_tif <- paste0(h5_name, ".tif")
    match_idx <- match(match_tif, patch_names)
    if (is.na(match_idx)) match_idx <- match(h5_name, patch_names)
    if (is.na(match_idx)) next

    patch_info <- dominant[match_idx, ]

    tryCatch({
      f <- hdf5r::H5File$new(h5_files[i], mode = "r")
      on.exit(f$close_all(), add = TRUE)

      # sen-2-data : (T, 10, 6, 6)
      if (f$exists("sen-2-data")) {
        s2 <- f[["sen-2-data"]]$read()
        # Spatial mean over 6x6 → (T, 10)
        n_t <- dim(s2)[1]
        n_b <- dim(s2)[2]

        # Moyenne spatiale
        s2_mean <- apply(s2, c(1, 2), mean, na.rm = TRUE)

        # Interpoler à 73 dates si nécessaire
        if (n_t != length(target_dates)) {
          orig_dates <- seq.Date(
            as.Date(TS_PARAMS$start_date),
            as.Date(TS_PARAMS$end_date),
            length.out = n_t
          )
          s2_interp <- matrix(NA, nrow = length(target_dates), ncol = n_b)
          for (b in 1:n_b) {
            s2_interp[, b] <- approx(
              as.numeric(orig_dates), s2_mean[, b],
              xout = as.numeric(target_dates),
              rule = 2
            )$y
          }
          s2_mean <- s2_interp
          n_t <- length(target_dates)
        }

        # Construire le data.frame (1 ligne par date)
        bands_df <- as.data.frame(s2_mean)
        if (ncol(bands_df) >= 10) {
          names(bands_df) <- S2_BAND_NAMES[1:ncol(bands_df)]
        }

        patch_df <- data.frame(
          plot_id      = patch_info$patch_id,
          species_code = patch_info$species_code,
          species_name = patch_info$species_name,
          date         = target_dates[1:n_t],
          stringsAsFactors = FALSE
        )
        patch_df <- cbind(patch_df, bands_df)

        # Calculer les indices spectraux
        if (all(c("B08", "B04", "B02") %in% names(patch_df))) {
          patch_df$NDVI   <- calc_ndvi(patch_df$B08, patch_df$B04)
          patch_df$EVI    <- calc_evi(patch_df$B08, patch_df$B04, patch_df$B02)
          patch_df$NDWI   <- (patch_df$B03 - patch_df$B08) / (patch_df$B03 + patch_df$B08 + 1e-10)
          patch_df$NBR    <- (patch_df$B08 - patch_df$B12) / (patch_df$B08 + patch_df$B12 + 1e-10)
          patch_df$CRI    <- (1 / (patch_df$B02 + 1e-10)) - (1 / (patch_df$B03 + 1e-10))
          patch_df$RENDVI <- (patch_df$B06 - patch_df$B05) / (patch_df$B06 + patch_df$B05 + 1e-10)
        }

        all_rows[[length(all_rows) + 1]] <- patch_df
      }
    }, error = function(e) {
      # Ignorer les fichiers illisibles
    })
  }

  cli::cli_progress_done()

  if (length(all_rows) == 0) {
    cli::cli_alert_warning("Aucun patch HDF5 lisible")
    return(NULL)
  }

  result <- do.call(rbind, all_rows)
  cli::cli_alert_success("{length(all_rows)} patchs lus, {nrow(result)} observations")
  result
}

# --- Génération synthétique basée sur les labels réels ---
.generate_ts_from_labels <- function(dominant) {
  target_dates <- seq.Date(
    as.Date(TS_PARAMS$start_date),
    as.Date(TS_PARAMS$end_date),
    by = TS_PARAMS$target_interval_days
  )
  n_dates  <- length(target_dates)
  n_patches <- nrow(dominant)
  total_rows <- n_patches * n_dates

  cli::cli_alert_info("G\u00e9n\u00e9ration pour {n_patches} patchs \u00d7 {n_dates} dates ({total_rows} lignes)")

  # Pré-allouer les vecteurs de métadonnées
  plot_ids <- rep(dominant$patch_id, each = n_dates)
  sp_codes <- rep(dominant$species_code, each = n_dates)
  sp_names <- rep(dominant$species_name, each = n_dates)
  dates_vec <- rep(target_dates, times = n_patches)

  # Hash du nom de genre → offset pour varier les profils entre genres
  # partageant le même species_code proxy
  genus_offset <- as.integer(
    vapply(dominant$species_name, function(g) sum(utf8ToInt(g)), numeric(1))
  )

  # Pré-allouer la matrice de bandes (beaucoup plus rapide que 50k data.frames)
  band_mat <- matrix(0, nrow = total_rows, ncol = 10)

  cli::cli_progress_bar(
    "G\u00e9n\u00e9ration synth\u00e9tique",
    total = n_patches,
    format = "{cli::pb_bar} {cli::pb_percent} | {cli::pb_current}/{cli::pb_total} patchs | ETA: {cli::pb_eta}"
  )

  for (i in seq_len(n_patches)) {
    if (i %% 100 == 0 || i == n_patches) cli::cli_progress_update(set = i)

    sp_code <- dominant$species_code[i]
    if (sp_code == 0 || is.na(sp_code)) sp_code <- 1L

    # Graine unique par genre + patch (même code proxy → profils distincts grâce à l'offset)
    set.seed(sp_code * 1000L + genus_offset[i] + i)

    # Simuler NDVI via le moteur existant
    ndvi <- simulate_species_ndvi(sp_code, target_dates)

    # Dériver les bandes spectrales
    nir   <- 0.3 + 0.4 * ndvi + rnorm(n_dates, 0, 0.02)
    red   <- nir * (1 - ndvi) / (1 + ndvi + 1e-6) + rnorm(n_dates, 0, 0.01)
    blue  <- red * runif(1, 0.7, 0.9) + rnorm(n_dates, 0, 0.01)
    green <- (red + nir) / 3 + rnorm(n_dates, 0, 0.01)
    re1   <- (red + nir) / 2 * runif(1, 0.85, 0.95) + rnorm(n_dates, 0, 0.01)
    re2   <- (re1 + nir) / 2 + rnorm(n_dates, 0, 0.01)
    re3   <- nir * runif(1, 0.9, 0.98) + rnorm(n_dates, 0, 0.01)
    nir2  <- nir * runif(1, 0.85, 0.95) + rnorm(n_dates, 0, 0.01)
    swir1 <- 0.2 - 0.1 * ndvi + rnorm(n_dates, 0, 0.02)
    swir2 <- swir1 * runif(1, 0.6, 0.8) + rnorm(n_dates, 0, 0.01)

    rows <- ((i - 1L) * n_dates + 1L):(i * n_dates)
    band_mat[rows,  1] <- pmax(0, blue)
    band_mat[rows,  2] <- pmax(0, green)
    band_mat[rows,  3] <- pmax(0, red)
    band_mat[rows,  4] <- pmax(0, re1)
    band_mat[rows,  5] <- pmax(0, re2)
    band_mat[rows,  6] <- pmax(0, re3)
    band_mat[rows,  7] <- pmax(0, nir)
    band_mat[rows,  8] <- pmax(0, nir2)
    band_mat[rows,  9] <- pmax(0, swir1)
    band_mat[rows, 10] <- pmax(0, swir2)
  }

  cli::cli_progress_done()

  # Construction d'un seul data.frame (évite do.call(rbind, 50k))
  cli::cli_alert_info("Construction du data.frame ({total_rows} lignes)...")

  result <- data.frame(
    plot_id      = plot_ids,
    species_code = sp_codes,
    species_name = sp_names,
    date         = dates_vec,
    B02 = band_mat[, 1], B03 = band_mat[, 2], B04 = band_mat[, 3],
    B05 = band_mat[, 4], B06 = band_mat[, 5], B07 = band_mat[, 6],
    B08 = band_mat[, 7], B8A = band_mat[, 8],
    B11 = band_mat[, 9], B12 = band_mat[, 10],
    stringsAsFactors = FALSE
  )

  # Indices spectraux — vectorisés sur tout le data.frame d'un coup
  cli::cli_alert_info("Calcul des indices spectraux...")
  result$NDVI   <- calc_ndvi(result$B08, result$B04)
  result$EVI    <- calc_evi(result$B08, result$B04, result$B02)
  result$NDWI   <- (result$B03 - result$B08) / (result$B03 + result$B08 + 1e-10)
  result$NBR    <- (result$B08 - result$B12) / (result$B08 + result$B12 + 1e-10)
  result$CRI    <- (1 / (result$B02 + 1e-10)) - (1 / (result$B03 + 1e-10))
  result$RENDVI <- (result$B06 - result$B05) / (result$B06 + result$B05 + 1e-10)

  cli::cli_alert_success("Dataset synth\u00e9tique : {nrow(result)} lignes ({n_patches} patchs)")
  result
}

# %||% est défini dans 08_download_satellite.R
