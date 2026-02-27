#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — PRÉDICTION SPATIALE SUR UNE ZONE D'INTÉRÊT
#
# Workflow :
#   aoi.gpkg → téléchargement S2 → séries temporelles pixel par pixel
#            → extraction de features → classification → carte des essences
#
# Usage :
#   source("R/07_predict_aoi.R")
#   carte <- predict_species_map("aoi.gpkg")
#
# En ligne de commande :
#   Rscript R/07_predict_aoi.R --aoi mon_aoi.gpkg --year 2023
# ==============================================================================

# Tous les modules chargés via le package

# ==============================================================================
# 1. TÉLÉCHARGEMENT DES SÉRIES TEMPORELLES SENTINEL-2 SUR L'AOI
# ==============================================================================

#' Recherche des images Sentinel-2 couvrant une AOI via Planetary Computer
#'
#' @param aoi sf object — polygone de la zone d'intérêt
#' @param year Année d'analyse
#' @param max_cloud Couverture nuageuse max par scène (%)
#' @param output_dir Répertoire de sortie pour les images
#' @return Liste avec les métadonnées des scènes
search_s2_for_aoi <- function(aoi, year = 2021, max_cloud = 30,
                               output_dir = file.path(RAW_DIR, "sentinel2")) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  bbox <- as.numeric(sf::st_bbox(aoi_wgs84))

  log_msg("Recherche Sentinel-2 L2A sur l'AOI (Planetary Computer)")
  log_msg("  Bbox : [{round(bbox[1],4)}, {round(bbox[2],4)}, {round(bbox[3],4)}, {round(bbox[4],4)}]")
  log_msg("  Période : {year}-01-01 → {year}-12-31")
  log_msg("  Nuages max : {max_cloud}%")

  start_date <- paste0(year, "-01-01")
  end_date   <- paste0(year, "-12-31")

  search_result <- search_sentinel2(aoi, start_date, end_date, max_cloud)
  if (is.null(search_result)) return(NULL)

  scenes_df <- search_result$scenes
  log_msg("  {nrow(scenes_df)} scènes trouvées", level = "success")
  log_msg("  Couverture : {min(scenes_df$date)} → {max(scenes_df$date)}")
  log_msg("  Nuages moyen : {round(mean(scenes_df$cloud_cover, na.rm = TRUE), 1)}%")

  # Retourner la liste pour compatibilité avec le reste du pipeline
  lapply(seq_len(nrow(scenes_df)), function(i) {
    list(
      id          = scenes_df$id[i],
      date        = scenes_df$date[i],
      cloud_cover = scenes_df$cloud_cover[i],
      platform    = scenes_df$platform[i]
    )
  })
}

#' Construction d'un cube raster spatio-temporel à partir d'images locales S2
#'
#' @param s2_dir Répertoire contenant les images Sentinel-2 (structure SAFE ou GeoTIFF)
#' @param aoi sf object — zone d'intérêt
#' @param bands Bandes à charger
#' @param year Année
#' @param resolution Résolution cible en mètres (10 ou 20)
#' @return SpatRaster empilé (bandes × dates)
build_s2_cube <- function(s2_dir, aoi, bands = S2_BAND_NAMES,
                           year = 2021, resolution = 10) {
  log_msg("Construction du cube Sentinel-2 depuis {s2_dir}")

  # Reprojection AOI en Lambert-93 si nécessaire
  aoi_proj <- sf::st_transform(aoi, 2154)
  aoi_vect <- terra::vect(aoi_proj)
  aoi_ext  <- terra::ext(aoi_vect)

  # Lister toutes les images
  all_files <- list.files(s2_dir, pattern = "\\.(tif|TIF|jp2|JP2)$",
                          recursive = TRUE, full.names = TRUE)

  if (length(all_files) == 0) {
    cli::cli_alert_danger("Aucune image trouvée dans {s2_dir}")
    return(NULL)
  }

  log_msg("  {length(all_files)} fichiers raster trouvés")

  # Organiser par date et bande
  file_info <- data.frame(
    path = all_files,
    filename = basename(all_files),
    stringsAsFactors = FALSE
  )

  # Extraire la date depuis le nom de fichier
  date_pattern <- "(\\d{4})(\\d{2})(\\d{2})"
  file_info$date_str <- stringr::str_extract(file_info$filename, date_pattern)
  file_info$date <- as.Date(file_info$date_str, format = "%Y%m%d")

  # Extraire la bande
  file_info$band <- stringr::str_extract(file_info$filename, "B\\d{2}|B8A|SCL")

  # Filtrer par année et bandes demandées
  file_info <- file_info[
    !is.na(file_info$date) &
    format(file_info$date, "%Y") == as.character(year) &
    file_info$band %in% c(bands, "SCL"),
  ]

  dates_available <- sort(unique(file_info$date))
  log_msg("  {length(dates_available)} dates disponibles en {year}")

  # Construire le cube : pour chaque date, empiler les bandes
  cube_list <- list()

  for (d in seq_along(dates_available)) {
    current_date <- dates_available[d]
    date_files <- file_info[file_info$date == current_date, ]

    band_stack <- list()
    for (band in bands) {
      bf <- date_files[date_files$band == band, ]
      if (nrow(bf) > 0) {
        r <- terra::rast(bf$path[1])

        # Reprojeter et rééchantillonner si nécessaire
        if (!terra::same.crs(r, terra::crs("EPSG:2154"))) {
          r <- terra::project(r, "EPSG:2154", method = "bilinear")
        }

        # Rééchantillonner à la résolution cible
        if (terra::res(r)[1] != resolution) {
          template <- terra::rast(aoi_ext, resolution = resolution, crs = "EPSG:2154")
          r <- terra::resample(r, template, method = "bilinear")
        }

        # Cropper sur l'AOI
        r <- terra::crop(r, aoi_ext)
        r <- terra::mask(r, aoi_vect)

        names(r) <- paste0(band, "_", format(current_date, "%Y%m%d"))
        band_stack[[band]] <- r
      }
    }

    if (length(band_stack) == length(bands)) {
      cube_list[[as.character(current_date)]] <- terra::rast(band_stack)
    }
  }

  log_msg("  Cube construit : {length(cube_list)} dates × {length(bands)} bandes",
          level = "success")

  cube_list
}

# ==============================================================================
# 2. CLASSIFICATION PIXEL PAR PIXEL SUR L'AOI
# ==============================================================================

#' Extraction des features pour chaque pixel du cube S2
#' (traitement par blocs pour gérer la mémoire)
#'
#' @param cube_list Liste de SpatRaster (1 par date, chaque couche = 1 bande)
#' @param dates Vecteur de dates correspondant au cube
#' @param block_size Nombre de lignes par bloc de traitement
#' @return SpatRaster multi-couches avec les features par pixel
extract_pixel_features <- function(cube_list, dates, block_size = 100) {
  log_msg("Extraction des features pixel par pixel...")

  n_dates  <- length(cube_list)
  template <- cube_list[[1]][[1]]  # Premier raster pour les dimensions
  n_cols   <- terra::ncol(template)
  n_rows   <- terra::nrow(template)
  n_pixels <- n_cols * n_rows

  log_msg("  Grille : {n_rows} × {n_cols} pixels ({n_pixels} pixels)")
  log_msg("  Dates : {n_dates}")

  # Dates cibles pour l'interpolation
  target_dates <- seq.Date(
    as.Date(paste0(format(dates[1], "%Y"), "-01-01")),
    as.Date(paste0(format(dates[1], "%Y"), "-12-31")),
    by = TS_PARAMS$target_interval_days
  )

  # Lire tout le cube en mémoire (par bande × date = matrice n_pixels × n_dates)
  bands <- S2_BAND_NAMES
  cube_arrays <- list()

  for (band in bands) {
    band_mat <- matrix(NA_real_, nrow = n_pixels, ncol = n_dates)
    for (d in seq_len(n_dates)) {
      layer_name <- grep(paste0("^", band, "_"), names(cube_list[[d]]), value = TRUE)
      if (length(layer_name) > 0) {
        vals <- terra::values(cube_list[[d]][[layer_name[1]]])
        band_mat[, d] <- as.numeric(vals)
      }
    }
    cube_arrays[[band]] <- band_mat
  }

  log_msg("  Cube chargé en mémoire")

  # Identifier les pixels valides (pas 100% NA)
  valid_mask <- rowSums(!is.na(cube_arrays[[bands[1]]])) >= 5
  valid_idx  <- which(valid_mask)
  n_valid    <- length(valid_idx)

  log_msg("  {n_valid} pixels valides sur {n_pixels} ({round(n_valid/n_pixels*100,1)}%)")

  if (n_valid == 0) {
    cli::cli_alert_danger("Aucun pixel valide dans l'AOI")
    return(NULL)
  }

  # --- Extraction des features pour chaque pixel valide ---
  log_msg("  Extraction des features...")
  pb <- cli::cli_progress_bar("Pixels", total = n_valid)

  feature_list <- vector("list", n_valid)

  for (i in seq_len(n_valid)) {
    px_idx <- valid_idx[i]

    # Récupérer les séries temporelles de toutes les bandes pour ce pixel
    bands_ts <- list()
    for (band in bands) {
      raw_vals <- cube_arrays[[band]][px_idx, ]

      # Interpolation aux dates cibles
      interp_vals <- interpolate_ts(dates, raw_vals, target_dates)

      # Lissage Savitzky-Golay
      if (!all(is.na(interp_vals))) {
        interp_vals <- smooth_savgol(interp_vals)
      }

      bands_ts[[band]] <- interp_vals
    }

    features <- c()

    # 1. Statistiques des bandes brutes
    for (band in bands) {
      stats <- calc_temporal_stats(bands_ts[[band]], prefix = band)
      features <- c(features, stats)
    }

    # 2. Indices spectraux
    if (all(c("B08", "B04", "B02", "B03", "B05", "B11", "B12") %in% names(bands_ts))) {
      indices <- calc_all_indices(bands_ts)

      for (idx_name in names(indices)) {
        stats <- calc_temporal_stats(indices[[idx_name]], prefix = idx_name)
        features <- c(features, stats)
      }

      # 3. Phénologie NDVI
      ndvi_smooth <- smooth_savgol(indices$NDVI)
      pheno <- extract_phenometrics(ndvi_smooth, target_dates)
      features <- c(features, pheno)

      # 4. Fourier NDVI
      fourier <- fourier_features(ndvi_smooth, n_harmonics = 3)
      names(fourier) <- paste0("NDVI_", names(fourier))
      features <- c(features, fourier)
    }

    feature_list[[i]] <- features
    cli::cli_progress_update(id = pb)
  }

  cli::cli_progress_done(id = pb)

  # Assembler en matrice
  feature_names <- names(feature_list[[1]])
  feature_mat <- matrix(NA_real_, nrow = n_valid, ncol = length(feature_names))
  colnames(feature_mat) <- feature_names
  for (i in seq_len(n_valid)) {
    feature_mat[i, ] <- feature_list[[i]]
  }

  log_msg("  Features extraites : {n_valid} pixels × {ncol(feature_mat)} features",
          level = "success")

  list(
    features   = feature_mat,
    valid_idx  = valid_idx,
    n_rows     = n_rows,
    n_cols     = n_cols,
    template   = template,
    feature_names = feature_names
  )
}

#' Classification des pixels et production de la carte des essences
#'
#' @param pixel_features Sortie de extract_pixel_features
#' @param model Modèle pré-entraîné (ranger)
#' @return Liste avec le raster classifié et les probabilités
classify_pixels <- function(pixel_features, model) {
  log_msg("Classification des pixels...")

  features_df <- as.data.frame(pixel_features$features)

  # Aligner les colonnes avec le modèle
  model_cols <- model$feature_cols
  missing_cols <- setdiff(model_cols, names(features_df))
  present_cols <- intersect(model_cols, names(features_df))

  if (length(missing_cols) > 0) {
    log_msg("  {length(missing_cols)} features du modèle absentes → remplies par 0",
            level = "warning")
    for (col in missing_cols) {
      features_df[[col]] <- 0
    }
  }

  log_msg("  {length(present_cols)}/{length(model_cols)} features alignées")

  # Remplacement des NA
  for (col in model_cols) {
    na_mask <- is.na(features_df[[col]])
    if (any(na_mask)) {
      features_df[[col]][na_mask] <- median(features_df[[col]], na.rm = TRUE)
    }
    # Si tout est NA, mettre 0
    if (all(is.na(features_df[[col]]))) {
      features_df[[col]] <- 0
    }
  }

  # Prédiction
  pred <- predict(model, data = features_df[, model_cols])

  # Classe prédite
  predicted_class <- apply(pred$predictions, 1, function(row) {
    which.max(row)
  })

  # Probabilité max (confiance)
  max_proba <- apply(pred$predictions, 1, max)

  log_msg("  Classification terminée", level = "success")

  list(
    class_idx    = predicted_class,
    class_names  = model$class_names,
    max_proba    = max_proba,
    all_probas   = pred$predictions,
    valid_idx    = pixel_features$valid_idx
  )
}

#' Reconstruction du raster classifié à partir des prédictions pixellaires
#'
#' @param predictions Sortie de classify_pixels
#' @param pixel_features Sortie de extract_pixel_features
#' @return Liste avec le raster d'espèces et le raster de confiance
build_species_raster <- function(predictions, pixel_features) {
  log_msg("Construction du raster des essences...")

  template <- pixel_features$template
  n_pixels <- pixel_features$n_rows * pixel_features$n_cols

  # Raster des classes (code numérique 1-20)
  class_vals <- rep(NA_real_, n_pixels)
  class_vals[predictions$valid_idx] <- predictions$class_idx

  r_class <- terra::rast(template)
  terra::values(r_class) <- class_vals
  names(r_class) <- "species_code"

  # Raster de confiance (probabilité max)
  proba_vals <- rep(NA_real_, n_pixels)
  proba_vals[predictions$valid_idx] <- predictions$max_proba

  r_proba <- terra::rast(template)
  terra::values(r_proba) <- proba_vals
  names(r_proba) <- "confidence"

  # Table d'attribution (code → nom d'espèce)
  levels_df <- data.frame(
    value  = seq_along(predictions$class_names),
    species = predictions$class_names
  )
  levels(r_class) <- levels_df

  log_msg("  Raster classifié construit", level = "success")

  list(species = r_class, confidence = r_proba, legend = levels_df)
}

# ==============================================================================
# 3. PIPELINE COMPLET : AOI → CARTE DES ESSENCES
# ==============================================================================

#' Pipeline complet de prédiction sur une zone d'intérêt
#'
#' @param aoi_path Chemin vers le fichier GeoPackage (ou Shapefile) de l'AOI
#' @param s2_dir Répertoire contenant les images Sentinel-2 L2A (NULL = auto-download)
#' @param s1_dir Répertoire contenant les images Sentinel-1 GRD (NULL = optionnel)
#' @param model_path Chemin vers le modèle pré-entraîné (.rds ou .pt)
#' @param year Année d'analyse
#' @param output_dir Répertoire de sortie
#' @param resolution Résolution cible en mètres
#' @param auto_download Télécharger automatiquement S2/S1 si pas de données locales
#' @param use_s1 Inclure les features Sentinel-1 dans la classification
#' @param use_pytorch Utiliser un modèle PyTorch (.pt) au lieu de ranger (.rds)
#'                    Nécessite l'environnement conda "treesat" (setup_python_env())
#' @param pytorch_model Architecture PyTorch : "tempcnn", "lstm", "transformer", "inception"
#' @return Liste avec les rasters et statistiques
predict_species_map <- function(aoi_path,
                                 s2_dir       = NULL,
                                 s1_dir       = NULL,
                                 model_path   = NULL,
                                 year         = 2021,
                                 output_dir   = OUTPUT_DIR,
                                 resolution   = 10,
                                 auto_download = FALSE,
                                 use_s1       = FALSE,
                                 use_pytorch  = FALSE,
                                 pytorch_model = "tempcnn") {

  cli::cli_h1("TreeSatAI — Détection d'essences sur zone d'intérêt")
  t_start <- Sys.time()

  # --- 1. Charger l'AOI ---
  cli::cli_h2("1. Chargement de la zone d'intérêt")
  aoi <- sf::st_read(aoi_path, quiet = TRUE)

  # Fusionner si multi-polygones
  if (nrow(aoi) > 1) {
    log_msg("  {nrow(aoi)} entités → fusion en un seul polygone")
    aoi <- sf::st_union(aoi) |> sf::st_as_sf()
  }

  aoi_area_ha <- as.numeric(sf::st_area(sf::st_transform(aoi, 2154))) / 10000
  log_msg("  Surface : {round(aoi_area_ha, 1)} ha")
  log_msg("  CRS : {sf::st_crs(aoi)$input}")

  # Estimation du nombre de pixels
  n_pixels_est <- round(aoi_area_ha * 10000 / resolution^2)
  log_msg("  Pixels estimés (~{resolution}m) : {format(n_pixels_est, big.mark = ' ')}")

  if (n_pixels_est > 5e6) {
    cli::cli_alert_warning("Zone très grande ({format(n_pixels_est, big.mark=' ')} pixels).")
    cli::cli_text("Considérez une résolution de 20m ou découpez la zone.")
  }

  # --- 2. Charger ou entraîner le modèle ---
  cli::cli_h2("2. Modèle de classification")

  py_model <- NULL  # Modèle PyTorch (si use_pytorch)

  if (use_pytorch) {
    # --- Mode PyTorch ---
    log_msg("  Mode : PyTorch ({pytorch_model})")

    if (!check_python_ready()) {
      cli::cli_alert_info("Initialisation de l'environnement Python...")
      setup_ok <- setup_python_env()
      if (!setup_ok) {
        cli::cli_alert_danger("Impossible d'initialiser Python. Repli sur ranger.")
        use_pytorch <- FALSE
      }
    }

    if (use_pytorch) {
      if (!is.null(model_path) && file.exists(model_path) && grepl("\\.pt$", model_path)) {
        py_model <- py_load_model(model_path)
      } else {
        # Chercher un modèle .pt existant
        default_pt <- file.path(MODELS_DIR, paste0("treesatai_", pytorch_model, "_best.pt"))
        if (file.exists(default_pt)) {
          py_model <- py_load_model(default_pt)
        } else {
          cli::cli_alert_warning("Aucun modèle PyTorch trouvé.")
          cli::cli_text("Entraînez d'abord : {.code py_train_model('data.csv', model_type = '{pytorch_model}')}")
          cli::cli_text("Ou depuis le terminal : {.code python python/train.py --data data.csv --model {pytorch_model}}")
          cli::cli_text("")
          cli::cli_alert_info("Repli sur Random Forest (ranger)")
          use_pytorch <- FALSE
        }
      }
    }
  }

  if (!use_pytorch) {
    # --- Mode ranger (R) ---
    # Chercher un modèle existant : chemin explicite, puis défaut
    default_rds <- file.path(MODELS_DIR, "treesatai_rf.rds")

    if (!is.null(model_path) && file.exists(model_path) && grepl("\\.rds$", model_path)) {
      model <- readRDS(model_path)
      log_msg("  Modèle ranger chargé depuis {model_path}", level = "success")
      log_msg("  Classes : {model$n_classes} espèces")

    } else if (file.exists(default_rds)) {
      model <- readRDS(default_rds)
      log_msg("  Modèle ranger chargé depuis {default_rds}", level = "success")
      log_msg("  Classes : {model$n_classes} espèces")

    } else {
      # --- Auto-entraînement sur les données synthétiques TreeSatAI ---
      cli::cli_alert_warning("Aucun modèle pré-entraîné trouvé.")
      cli::cli_alert_info("Auto-entraînement d'un modèle Random Forest sur les données TreeSatAI synthétiques...")
      cli::cli_text("")

      # 1. Générer le dataset synthétique (20 espèces × 50 échantillons)
      ts_long <- generate_synthetic_dataset(n_samples_per_species = 50, year = year)

      # 2. Construire la matrice de features
      feature_matrix <- build_feature_matrix(ts_long)

      # 3. Split train/test
      split <- split_train_test(feature_matrix, target_col = "species_name")
      train_data <- split$train
      test_data  <- split$test

      # 4. Sélectionner les features
      feature_cols <- select_features(feature_matrix)

      # 5. Entraîner le Random Forest
      model <- train_random_forest(train_data, feature_cols)

      # 6. Évaluer
      rf_preds <- predict_rf(model, test_data)
      rf_eval <- evaluate_classification(
        y_true      = test_data$species_name,
        y_pred      = rf_preds$predicted_class,
        class_names = SPECIES$french
      )

      # 7. Sauvegarder pour les prochaines utilisations
      save_model(model, rf_eval, model_name = "treesatai_rf")

      cli::cli_alert_success(
        "Modèle auto-entraîné — OA : {round(rf_eval$overall_accuracy * 100, 1)}%, Kappa : {round(rf_eval$kappa, 3)}"
      )
      cli::cli_alert_info(
        "Modèle sauvegardé dans {default_rds} — il sera réutilisé aux prochains appels."
      )
      cli::cli_text("")
      cli::cli_alert_warning(
        "Ce modèle est basé sur des données synthétiques. Pour de meilleurs résultats, entraînez sur des données réelles :"
      )
      cli::cli_text("  {.code Rscript R/06_pipeline.R --data /chemin/vers/parcelles}")
    }
  }

  # --- 3. Construire le cube Sentinel-2 ---
  cli::cli_h2("3. Données Sentinel-2")

  if (!is.null(s2_dir) && dir.exists(s2_dir)) {
    # Mode données locales
    cube_list <- build_s2_cube(s2_dir, aoi, bands = S2_BAND_NAMES,
                                year = year, resolution = resolution)

    if (is.null(cube_list) || length(cube_list) == 0) {
      cli::cli_alert_danger("Impossible de construire le cube S2")
      return(NULL)
    }

    dates <- as.Date(names(cube_list))

  } else if (auto_download) {
    # --- Mode téléchargement automatique via Planetary Computer (gratuit) ---
    cli::cli_h3("Téléchargement automatique (Planetary Computer — sans authentification)")

    sat_data <- download_satellite_data(
      aoi_path    = aoi_path,
      year        = year,
      download_s2 = TRUE,
      download_s1 = use_s1,
      max_cloud   = TS_PARAMS$max_cloud_cover,
      output_dir  = RAW_DIR
    )

    if (!is.null(sat_data$s2_dir)) {
      s2_bands_dir <- file.path(sat_data$s2_dir, "bands")
      if (!dir.exists(s2_bands_dir)) s2_bands_dir <- sat_data$s2_dir

      cube_list <- build_s2_cube(s2_bands_dir, aoi, bands = S2_BAND_NAMES,
                                  year = year, resolution = resolution)
      dates <- as.Date(names(cube_list))

      # Sentinel-1 RTC : déjà en dB, pas besoin de preprocess_s1()
      if (use_s1 && !is.null(sat_data$s1_dir)) {
        s1_bands_dir <- file.path(sat_data$s1_dir, "bands")
        if (!dir.exists(s1_bands_dir)) s1_bands_dir <- sat_data$s1_dir
        s1_cube <- build_s1_cube(s1_bands_dir, aoi, year, resolution)
      }
    } else {
      cli::cli_alert_danger("Téléchargement S2 échoué")
      cli::cli_text("Vérifiez votre connexion internet.")
      return(NULL)
    }

  } else {
    cli::cli_alert_danger("Pas de données satellite.")
    cli::cli_text("")
    cli::cli_text("Deux options :")
    cli::cli_text("")
    cli::cli_h3("Option A : Données locales")
    cli::cli_text('  predict_species_map("{aoi_path}", s2_dir = "/chemin/vers/S2/")')
    cli::cli_text("")
    cli::cli_h3("Option B : Téléchargement automatique (gratuit, sans compte)")
    cli::cli_text('  predict_species_map("{aoi_path}", auto_download = TRUE, year = {year})')
    cli::cli_text('  # Avec Sentinel-1 (radar) en plus :')
    cli::cli_text('  predict_species_map("{aoi_path}", auto_download = TRUE, use_s1 = TRUE)')
    stop("Fournissez s2_dir ou utilisez auto_download = TRUE")
  }

  # --- 4. Extraction des features ---
  cli::cli_h2("4. Extraction des features pixellaires")
  pixel_features <- extract_pixel_features(cube_list, dates)

  if (is.null(pixel_features)) {
    cli::cli_alert_danger("Échec de l'extraction des features")
    return(NULL)
  }

  # Intégration des features Sentinel-1 si disponibles
  if (use_s1 && exists("s1_cube") && !is.null(s1_cube)) {
    cli::cli_h3("4b. Features Sentinel-1 (radar)")
    s1_features <- extract_s1_pixel_features(s1_cube, pixel_features)
    if (!is.null(s1_features)) {
      # Fusionner avec les features S2
      pixel_features$features <- cbind(pixel_features$features, s1_features)
      pixel_features$feature_names <- colnames(pixel_features$features)
      log_msg("  Features S1 ajoutées : {ncol(s1_features)} colonnes radar")
    }
  }

  # --- 5. Classification ---
  cli::cli_h2("5. Classification des pixels")

  if (use_pytorch && !is.null(py_model)) {
    # --- Classification PyTorch ---
    predictions <- classify_pixels_pytorch(
      pixel_features, model_path = NULL,
      cube_list = cube_list, batch_size = 512L
    )
    # classify_pixels_pytorch charge le modèle — on le passe directement ici
    # pour éviter de le recharger :
    pixel_data <- cube_to_pytorch_array(
      cube_list, pixel_features$valid_idx, bands = S2_BAND_NAMES
    )
    py_results <- py_predict_pixels(py_model, pixel_data, batch_size = 512L)
    predictions <- list(
      class_idx    = py_results$predicted_class,
      class_names  = py_model$class_names,
      max_proba    = py_results$max_proba,
      all_probas   = py_results$probabilities,
      valid_idx    = pixel_features$valid_idx
    )
  } else {
    # --- Classification ranger (R) ---
    predictions <- classify_pixels(pixel_features, model)
  }

  # --- 6. Construction du raster ---
  cli::cli_h2("6. Production de la carte")
  rasters <- build_species_raster(predictions, pixel_features)

  # --- 7. Sauvegarde des résultats ---
  cli::cli_h2("7. Sauvegarde")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Raster GeoTIFF
  tif_path <- file.path(output_dir, "carte_essences.tif")
  terra::writeRaster(rasters$species, tif_path, overwrite = TRUE,
                     datatype = "INT1U")
  log_msg("  Raster espèces  : {tif_path}", level = "success")

  # Raster confiance
  proba_path <- file.path(output_dir, "carte_confiance.tif")
  terra::writeRaster(rasters$confidence, proba_path, overwrite = TRUE)
  log_msg("  Raster confiance : {proba_path}", level = "success")

  # Vectoriser le raster (polygones par espèce)
  log_msg("  Vectorisation...")
  species_poly <- terra::as.polygons(rasters$species, dissolve = TRUE)
  species_sf   <- sf::st_as_sf(species_poly)

  # Ajouter les noms d'espèces
  if ("species_code" %in% names(species_sf)) {
    species_sf <- species_sf |>
      dplyr::left_join(
        data.frame(species_code = rasters$legend$value,
                   species_name = rasters$legend$species),
        by = "species_code"
      )
  }

  gpkg_path <- file.path(output_dir, "carte_essences.gpkg")
  sf::st_write(species_sf, gpkg_path, delete_dsn = TRUE, quiet = TRUE)
  log_msg("  Vecteur espèces  : {gpkg_path}", level = "success")

  # Légende CSV
  legend_path <- file.path(output_dir, "legende_especes.csv")
  readr::write_csv(rasters$legend, legend_path)

  # Statistiques
  stats <- compute_map_statistics(predictions, rasters)
  stats_path <- file.path(output_dir, "statistiques_essences.csv")
  readr::write_csv(stats, stats_path)
  log_msg("  Statistiques     : {stats_path}", level = "success")

  # --- 8. Résumé ---
  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  cli::cli_h1("Résultat")
  cli::cli_text("")
  cli::cli_alert_success("Carte des essences produite en {round(t_elapsed, 1)} minutes")
  cli::cli_text("")

  print_species_summary(stats)

  cli::cli_text("")
  cli::cli_text("Fichiers de sortie :")
  cli::cli_ul()
  cli::cli_li("{tif_path}   — raster classifié (GeoTIFF)")
  cli::cli_li("{proba_path} — confiance de prédiction")
  cli::cli_li("{gpkg_path}  — polygones par espèce (GeoPackage)")
  cli::cli_li("{stats_path} — statistiques par espèce")
  cli::cli_end()

  invisible(list(
    species_raster    = rasters$species,
    confidence_raster = rasters$confidence,
    species_vector    = species_sf,
    statistics        = stats,
    model             = model,
    legend            = rasters$legend
  ))
}

# ==============================================================================
# 3b. EXTRACTION DES FEATURES SENTINEL-1 PIXEL PAR PIXEL
# ==============================================================================

#' Extraction des features radar S1 pour les mêmes pixels que S2
#' @param s1_cube Liste de SpatRaster S1 (VV, VH, ratio par date)
#' @param pixel_features Résultat de extract_pixel_features (pour valid_idx, template)
#' @return Matrice de features S1 (n_valid × n_s1_features)
extract_s1_pixel_features <- function(s1_cube, pixel_features) {
  log_msg("Extraction des features Sentinel-1...")

  if (is.null(s1_cube) || length(s1_cube) == 0) {
    log_msg("  Pas de données S1 disponibles", level = "warning")
    return(NULL)
  }

  n_dates_s1 <- length(s1_cube)
  valid_idx  <- pixel_features$valid_idx
  n_valid    <- length(valid_idx)
  n_pixels   <- pixel_features$n_rows * pixel_features$n_cols

  log_msg("  {n_dates_s1} dates S1 × {n_valid} pixels valides")

  # Charger le cube S1 en mémoire
  s1_pols <- c("VV", "VH")
  s1_arrays <- list()

  for (pol in s1_pols) {
    pol_mat <- matrix(NA_real_, nrow = n_pixels, ncol = n_dates_s1)
    for (d in seq_len(n_dates_s1)) {
      layer_name <- grep(paste0("^", pol, "_"), names(s1_cube[[d]]), value = TRUE)
      if (length(layer_name) > 0) {
        vals <- terra::values(s1_cube[[d]][[layer_name[1]]])
        pol_mat[, d] <- as.numeric(vals)
      }
    }
    s1_arrays[[pol]] <- pol_mat
  }

  # Extraire les features S1 par pixel
  feature_list <- vector("list", n_valid)

  for (i in seq_len(n_valid)) {
    px <- valid_idx[i]
    vv_ts <- s1_arrays$VV[px, ]
    vh_ts <- s1_arrays$VH[px, ]

    feature_list[[i]] <- calc_s1_temporal_features(vv_ts, vh_ts)
  }

  # Assembler en matrice
  feature_names <- names(feature_list[[1]])
  feature_mat <- matrix(NA_real_, nrow = n_valid, ncol = length(feature_names))
  colnames(feature_mat) <- feature_names
  for (i in seq_len(n_valid)) {
    feature_mat[i, ] <- feature_list[[i]]
  }

  log_msg("  Features S1 : {ncol(feature_mat)} colonnes", level = "success")
  feature_mat
}

# ==============================================================================
# 4. FONCTIONS AUXILIAIRES
# ==============================================================================

#' Calcul des statistiques de la carte classifiée
compute_map_statistics <- function(predictions, rasters) {
  class_names <- predictions$class_names
  class_idx   <- predictions$class_idx

  stats <- data.frame(
    code         = seq_along(class_names),
    espece       = class_names,
    n_pixels     = as.integer(table(factor(class_idx, levels = seq_along(class_names)))),
    stringsAsFactors = FALSE
  )

  total_pixels <- sum(stats$n_pixels)
  stats$pct <- round(stats$n_pixels / total_pixels * 100, 2)

  # Surface estimée (en hectares)
  res <- terra::res(rasters$species)[1]
  stats$surface_ha <- round(stats$n_pixels * res^2 / 10000, 2)

  # Confiance moyenne par espèce
  proba_by_class <- tapply(predictions$max_proba, class_idx, mean, na.rm = TRUE)
  stats$confiance_moy <- round(as.numeric(proba_by_class[as.character(stats$code)]) * 100, 1)
  stats$confiance_moy[is.na(stats$confiance_moy)] <- 0

  # Trier par surface décroissante
  stats <- stats[order(-stats$n_pixels), ]

  # Ajouter le type
  stats <- stats |>
    dplyr::left_join(SPECIES[, c("french", "type", "phenologie")],
                     by = c("espece" = "french"))

  stats
}

#' Affichage du résumé des espèces détectées
print_species_summary <- function(stats) {
  # Filtrer les espèces détectées
  detected <- stats[stats$n_pixels > 0, ]

  cli::cli_h2("Espèces détectées : {nrow(detected)} / {nrow(stats)}")

  if (nrow(detected) == 0) {
    cli::cli_alert_warning("Aucune espèce détectée")
    return(invisible(NULL))
  }

  cli::cli_text("")

  # Top espèces
  top_n <- min(10, nrow(detected))
  for (i in 1:top_n) {
    sp <- detected[i, ]
    bar_len <- round(sp$pct / max(detected$pct) * 30)
    bar <- paste0(rep("\u2588", bar_len), collapse = "")
    type_icon <- if (!is.na(sp$type) && sp$type == "feuillu") "\U0001F333" else "\U0001F332"

    cli::cli_text("  {type_icon} {sprintf('%-25s', sp$espece)} {bar} {sp$pct}% ({sp$surface_ha} ha)")
  }

  if (nrow(detected) > top_n) {
    cli::cli_text("  ... et {nrow(detected) - top_n} autres espèces")
  }

  # Résumé par type
  cli::cli_text("")
  by_type <- aggregate(surface_ha ~ type, data = detected, sum, na.rm = TRUE)
  for (i in seq_len(nrow(by_type))) {
    cli::cli_text("  {by_type$type[i]} : {by_type$surface_ha[i]} ha")
  }

  invisible(detected)
}

# ==============================================================================
# 5. EXÉCUTION EN LIGNE DE COMMANDE
# Déplacé dans inst/scripts/07_predict_cli.R pour compatibilité package
# Usage CLI : Rscript inst/scripts/07_predict_cli.R --aoi aoi.gpkg
# ==============================================================================

