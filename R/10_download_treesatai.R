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
        cli::cli_alert_info("  D\u00e9j\u00e0 pr\u00e9sent : {file_name}")
        next
      }

      # URL de téléchargement direct
      dl_url <- paste0(HF_RESOLVE_BASE, "/", HF_REPO_ID, "/resolve/main/", file_path)

      size_mb <- round(file_size / 1024 / 1024, 1)
      cli::cli_alert_info("  T\u00e9l\u00e9chargement : {file_name} ({size_mb} Mo)")

      tryCatch({
        download.file(dl_url, dest_file, mode = "wb", quiet = TRUE)
        cli::cli_alert_success("  {file_name} OK")
      }, error = function(e) {
        cli::cli_alert_danger("  \u00c9chec : {e$message}")
      })
    }
  }

  # --- Post-traitement : extraction des zips si présents ---
  zip_files <- list.files(dest_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
  if (length(zip_files) > 0) {
    cli::cli_h2("Extraction des archives")
    for (zf in zip_files) {
      cli::cli_alert_info("Extraction : {basename(zf)}")
      tryCatch({
        unzip(zf, exdir = dirname(zf))
        cli::cli_alert_success("  Extrait dans {dirname(zf)}")
      }, error = function(e) {
        cli::cli_alert_warning("  \u00c9chec extraction : {e$message}")
      })
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
