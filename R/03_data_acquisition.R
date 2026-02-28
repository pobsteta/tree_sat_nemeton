#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Acquisition et préparation des données
# Téléchargement et structuration des séries temporelles Sentinel-2
# ==============================================================================

# Configuration et utilitaires chargés via le package

# ==============================================================================
# OPTION 1 : Chargement du dataset TreeSatAI-Time-Series pré-constitué
# ==============================================================================

#' Téléchargement du dataset TreeSatAI-Time-Series depuis Zenodo/IGN
#' @param dest_dir Répertoire de destination
#' @param source URL ou chemin local du dataset
#' @return Chemin vers les données téléchargées
download_treesatai_ts <- function(dest_dir = RAW_DIR, source = NULL) {
  log_msg("Acquisition du dataset TreeSatAI-Time-Series (IGNF)")

  # Le dataset est distribué par l'IGN France
  # URL type Zenodo (à adapter selon la publication)
  if (is.null(source)) {
    cli::cli_alert_info(paste(
      "Le dataset TreeSatAI-Time-Series n'est pas directement téléchargeable.",
      "Veuillez le placer manuellement dans : {dest_dir}",
      "",
      "Structure attendue :",
      "  {dest_dir}/",
      "    ├── plots/          # Parcelles de référence (GeoPackage/Shapefile)",
      "    ├── timeseries/     # Séries temporelles extraites (.csv / .parquet)",
      "    └── metadata.json   # Métadonnées du dataset",
      sep = "\n"
    ))
    return(invisible(dest_dir))
  }

  # Si une URL est fournie, téléchargement
  dest_file <- file.path(dest_dir, "treesatai_ts.zip")
  if (!file.exists(dest_file)) {
    log_msg("Téléchargement depuis {source}...")
    download.file(source, dest_file, mode = "wb")
    log_msg("Extraction de l'archive...", level = "info")
    unzip(dest_file, exdir = dest_dir)
    log_msg("Dataset extrait dans {dest_dir}", level = "success")
  } else {
    log_msg("Archive déjà présente : {dest_file}", level = "info")
  }

  invisible(dest_dir)
}

# ==============================================================================
# OPTION 2 : Construction du dataset à partir de l'IFN et Sentinel-2
# ==============================================================================

#' Chargement des placettes IFN (Inventaire Forestier National)
#' @param ifn_path Chemin vers le fichier des placettes IFN
#' @param species_filter Vecteur de codes espèces à conserver (optionnel)
#' @return sf object avec les placettes filtrées
load_ifn_plots <- function(ifn_path, species_filter = NULL) {
  log_msg("Chargement des placettes IFN...")

  plots <- sf::st_read(ifn_path, quiet = TRUE)
  log_msg("  {nrow(plots)} placettes chargées", level = "info")

  # Vérification des colonnes requises
  required_cols <- c("idp", "essence", "geometry")
  missing_cols <- setdiff(required_cols, tolower(names(plots)))
  if (length(missing_cols) > 0) {
    cli::cli_alert_warning("Colonnes manquantes : {paste(missing_cols, collapse = ', ')}")
    cli::cli_text("Colonnes disponibles : {paste(names(plots), collapse = ', ')}")
  }

  # Normaliser les noms de colonnes
  names(plots) <- tolower(names(plots))

  # Filtrage par espèce si demandé
  if (!is.null(species_filter)) {
    plots <- plots[plots$essence %in% species_filter, ]
    log_msg("  {nrow(plots)} placettes après filtrage ({length(species_filter)} espèces)",
            level = "info")
  }

  # Assurer la projection en WGS84 / EPSG:4326 ou Lambert-93 / EPSG:2154
  if (sf::st_crs(plots)$epsg != 2154) {
    plots <- sf::st_transform(plots, 2154)
    log_msg("  Reprojection en Lambert-93 (EPSG:2154)", level = "info")
  }

  plots
}

#' Extraction des séries temporelles Sentinel-2 pour un ensemble de parcelles
#' @param plots sf object — points ou polygones des parcelles
#' @param s2_dir Répertoire contenant les images Sentinel-2 L2A
#' @param bands Bandes spectrales à extraire
#' @param buffer_m Rayon du buffer autour des points (en mètres)
#' @return data.frame avec colonnes : plot_id, date, band, value
extract_s2_timeseries <- function(plots, s2_dir, bands = S2_BAND_NAMES,
                                   buffer_m = 15) {
  log_msg("Extraction des séries temporelles Sentinel-2")
  log_msg("  {nrow(plots)} parcelles × {length(bands)} bandes", level = "info")

  # Créer un buffer si les parcelles sont des points
  if (all(sf::st_geometry_type(plots) == "POINT")) {
    plots_buf <- sf::st_buffer(plots, buffer_m)
    log_msg("  Buffer de {buffer_m}m appliqué aux points", level = "info")
  } else {
    plots_buf <- plots
  }

  # Lister les images par bande
  all_results <- list()

  for (band in bands) {
    log_msg("  Traitement bande {band}...")

    files_df <- list_s2_files(s2_dir, band)
    if (nrow(files_df) == 0) {
      cli::cli_alert_warning("Aucune image trouvée pour {band}")
      next
    }

    # Filtrage temporel
    files_df <- files_df[
      files_df$date >= as.Date(TS_PARAMS$start_date) &
      files_df$date <= as.Date(TS_PARAMS$end_date), ]

    log_msg("    {nrow(files_df)} images entre {TS_PARAMS$start_date} et {TS_PARAMS$end_date}",
            level = "info")

    # Extraction par image
    band_results <- lapply(seq_len(nrow(files_df)), function(i) {
      r <- terra::rast(files_df$path[i])

      # Extraction (moyenne sur le buffer/polygone)
      vals <- terra::extract(r, terra::vect(plots_buf), fun = mean, na.rm = TRUE)

      data.frame(
        plot_id = plots_buf$idp,
        date    = files_df$date[i],
        band    = band,
        value   = vals[, 2],  # 2ème colonne = valeurs
        stringsAsFactors = FALSE
      )
    })

    all_results[[band]] <- do.call(rbind, band_results)
  }

  result <- do.call(rbind, all_results)
  log_msg("Extraction terminée : {nrow(result)} observations", level = "success")

  result
}

#' Extraction en parallèle avec future
#' @param plots sf object
#' @param s2_dir Répertoire Sentinel-2
#' @param n_workers Nombre de workers
#' @return data.frame de séries temporelles
extract_s2_timeseries_parallel <- function(plots, s2_dir,
                                            bands = S2_BAND_NAMES,
                                            n_workers = 4) {
  log_msg("Extraction parallèle ({n_workers} workers)")

  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential))

  results <- future.apply::future_lapply(bands, function(band) {
    files_df <- list_s2_files(s2_dir, band)
    files_df <- files_df[
      files_df$date >= as.Date(TS_PARAMS$start_date) &
      files_df$date <= as.Date(TS_PARAMS$end_date), ]

    if (nrow(files_df) == 0) return(NULL)

    # Buffer si points
    plots_local <- if (all(sf::st_geometry_type(plots) == "POINT")) {
      sf::st_buffer(plots, 15)
    } else {
      plots
    }

    band_vals <- lapply(seq_len(nrow(files_df)), function(i) {
      r <- terra::rast(files_df$path[i])
      vals <- terra::extract(r, terra::vect(plots_local), fun = mean, na.rm = TRUE)
      data.frame(
        plot_id = plots_local$idp,
        date    = files_df$date[i],
        band    = band,
        value   = vals[, 2],
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, band_vals)
  }, future.seed = TRUE)

  result <- do.call(rbind, results)
  log_msg("Extraction parallèle terminée : {nrow(result)} observations", level = "success")
  result
}

# ==============================================================================
# OPTION 3 : Utilisation de l'API STAC pour Sentinel-2 (Copernicus Data Space)
# ==============================================================================

#' Recherche d'images Sentinel-2 via STAC API (Copernicus Data Space)
#' @param bbox Bounding box c(xmin, ymin, xmax, ymax) en EPSG:4326
#' @param start_date Date de début "YYYY-MM-DD"
#' @param end_date Date de fin "YYYY-MM-DD"
#' @param max_cloud Couverture nuageuse max (%)
#' @return data.frame des scènes trouvées
search_s2_stac <- function(bbox, start_date, end_date, max_cloud = 30) {
  log_msg("Recherche STAC Sentinel-2 L2A...")

  stac_url <- "https://catalogue.dataspace.copernicus.eu/stac"

  # Construction de la requête STAC
  body <- list(
    collections = list("sentinel-2-l2a"),
    bbox = bbox,
    datetime = paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z"),
    limit = 500,
    query = list(
      `eo:cloud_cover` = list(lte = max_cloud)
    )
  )

  resp <- httr2::request(paste0(stac_url, "/search")) |>
    httr2::req_body_json(body) |>
    httr2::req_perform()

  items <- httr2::resp_body_json(resp)
  n_items <- length(items$features)
  log_msg("  {n_items} scènes Sentinel-2 L2A trouvées", level = "success")

  if (n_items == 0) return(data.frame())

  # Extraction des métadonnées
  scenes <- lapply(items$features, function(feat) {
    data.frame(
      id         = feat$id,
      datetime   = feat$properties$datetime,
      cloud_cover = feat$properties$`eo:cloud_cover`,
      platform   = feat$properties$platform,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, scenes)
}

# ==============================================================================
# Préparation du jeu de données final
# ==============================================================================

#' Pivot des séries temporelles en format large (1 ligne = 1 parcelle)
#' @param ts_long data.frame long (plot_id, date, band, value)
#' @return data.frame large avec interpolation temporelle
prepare_ts_wide <- function(ts_long) {
  log_msg("Préparation du jeu de données en format large...")

  # Dates cibles pour l'interpolation
  target_dates <- seq.Date(
    as.Date(TS_PARAMS$start_date),
    as.Date(TS_PARAMS$end_date),
    by = TS_PARAMS$target_interval_days
  )

  # Liste des parcelles et bandes
  plot_ids <- unique(ts_long$plot_id)
  bands    <- unique(ts_long$band)

  log_msg("  {length(plot_ids)} parcelles × {length(bands)} bandes × {length(target_dates)} dates",
          level = "info")

  # Interpolation par parcelle et bande
  results <- lapply(plot_ids, function(pid) {
    plot_data <- ts_long[ts_long$plot_id == pid, ]

    row <- c(plot_id = pid)

    for (band in bands) {
      band_data <- plot_data[plot_data$band == band, ]
      band_data <- band_data[order(band_data$date), ]

      # Interpolation
      interp_values <- interpolate_ts(
        dates        = band_data$date,
        values       = band_data$value,
        target_dates = target_dates
      )

      # Lissage Savitzky-Golay
      smooth_values <- smooth_savgol(interp_values)

      # Nommage : band_DOY
      doy_labels <- format(target_dates, "%j")
      names(smooth_values) <- paste0(band, "_d", doy_labels)

      row <- c(row, smooth_values)
    }

    row
  })

  df <- as.data.frame(do.call(rbind, results), stringsAsFactors = FALSE)

  # Conversion en numérique (sauf plot_id)
  numeric_cols <- setdiff(names(df), "plot_id")
  df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

  log_msg("Dataset préparé : {nrow(df)} parcelles × {ncol(df)} features",
          level = "success")

  df
}

#' Ajout des indices spectraux et métriques phénologiques
#' @param ts_wide data.frame large (sortie de prepare_ts_wide)
#' @return data.frame enrichi
enrich_features <- function(ts_wide) {
  log_msg("Enrichissement avec indices spectraux et métriques phénologiques...")

  target_dates <- seq.Date(
    as.Date(TS_PARAMS$start_date),
    as.Date(TS_PARAMS$end_date),
    by = TS_PARAMS$target_interval_days
  )

  n_dates <- length(target_dates)
  doy_labels <- format(target_dates, "%j")

  enriched_rows <- lapply(seq_len(nrow(ts_wide)), function(i) {
    row <- ts_wide[i, ]

    # Reconstruire les bandes temporelles
    bands_ts <- list()
    for (band in S2_BAND_NAMES) {
      col_names <- paste0(band, "_d", doy_labels)
      existing <- col_names[col_names %in% names(row)]
      if (length(existing) > 0) {
        bands_ts[[band]] <- as.numeric(row[existing])
      }
    }

    extra_features <- c()

    # Indices spectraux (séries temporelles complètes)
    if (all(c("B08", "B04", "B02", "B03", "B05", "B11", "B12") %in% names(bands_ts))) {
      indices_df <- calc_all_indices(bands_ts)

      # Ajouter chaque indice comme série temporelle
      for (idx_name in names(indices_df)) {
        idx_vals <- indices_df[[idx_name]]
        idx_names <- paste0(idx_name, "_d", doy_labels[seq_along(idx_vals)])
        extra_features <- c(extra_features, setNames(idx_vals, idx_names))

        # Statistiques temporelles par indice
        stats <- calc_temporal_stats(idx_vals, prefix = idx_name)
        extra_features <- c(extra_features, stats)
      }

      # Métriques phénologiques (sur NDVI et EVI)
      ndvi_ts <- indices_df$NDVI
      evi_ts  <- indices_df$EVI

      ndvi_smooth <- smooth_savgol(ndvi_ts)
      evi_smooth  <- smooth_savgol(evi_ts)

      pheno_ndvi <- extract_phenometrics(ndvi_smooth, target_dates)
      pheno_evi  <- extract_phenometrics(evi_smooth, target_dates)
      names(pheno_evi) <- gsub("^pheno_", "EVI_pheno_", names(pheno_evi))

      extra_features <- c(extra_features, pheno_ndvi, pheno_evi)

      # Features Fourier
      fourier_ndvi <- fourier_features(ndvi_smooth, n_harmonics = 3)
      names(fourier_ndvi) <- paste0("NDVI_", names(fourier_ndvi))
      extra_features <- c(extra_features, fourier_ndvi)

      # Classification phénologique automatique
      pheno_type <- classify_phenotype(ndvi_smooth, target_dates)
      extra_features <- c(extra_features, pheno_class = match(pheno_type,
        c("evergreen", "deciduous", "semi_deciduous")))
    }

    # Statistiques sur les bandes brutes
    for (band in names(bands_ts)) {
      stats <- calc_temporal_stats(bands_ts[[band]], prefix = band)
      extra_features <- c(extra_features, stats)
    }

    extra_features
  })

  # Combiner
  extra_df <- as.data.frame(do.call(rbind, enriched_rows))
  result <- cbind(ts_wide, extra_df)

  log_msg("Features enrichies : {ncol(result)} colonnes au total", level = "success")
  result
}

# ==============================================================================
# Génération de données synthétiques (pour démonstration / test)
# ==============================================================================

#' Génération d'un profil NDVI typique par espèce
#' @param species_code Code espèce (1-20)
#' @param dates Vecteur de dates
#' @param noise_sd Écart-type du bruit gaussien
#' @return Vecteur NDVI simulé
simulate_species_ndvi <- function(species_code, dates, noise_sd = 0.03) {
  doy <- as.numeric(format(dates, "%j"))
  sp  <- SPECIES[species_code, ]

  set.seed(species_code * 100 + as.numeric(dates[1]))

  if (sp$phenologie == "sempervirent") {
    # Profil relativement plat, NDVI élevé toute l'année
    base_ndvi <- switch(sp$french,
      "Chêne vert"      = 0.65,
      "Épicéa commun"   = 0.55,
      "Sapin pectiné"    = 0.58,
      "Douglas"          = 0.60,
      "Pin sylvestre"    = 0.50,
      "Pin maritime"     = 0.52,
      "Pin noir"         = 0.48,
      "Pin d'Alep"       = 0.45,
      0.55  # défaut
    )
    # Légère variation saisonnière
    amplitude <- runif(1, 0.05, 0.12)
    phase_shift <- runif(1, -10, 10)
    ndvi <- base_ndvi + amplitude * sin(2 * pi * (doy + phase_shift) / 365)

  } else {
    # Profil caduc : forte amplitude saisonnière
    params <- switch(sp$french,
      "Chêne pédonculé"  = list(min = 0.20, max = 0.82, sos = 100, peak = 170, eos = 300),
      "Chêne sessile"    = list(min = 0.22, max = 0.80, sos = 105, peak = 175, eos = 295),
      "Chêne pubescent"  = list(min = 0.18, max = 0.78, sos = 95,  peak = 165, eos = 305),
      "Hêtre"            = list(min = 0.15, max = 0.85, sos = 110, peak = 180, eos = 290),
      "Châtaignier"      = list(min = 0.18, max = 0.83, sos = 105, peak = 175, eos = 295),
      "Charme"           = list(min = 0.17, max = 0.80, sos = 100, peak = 170, eos = 300),
      "Bouleau verruqueux" = list(min = 0.12, max = 0.78, sos = 90, peak = 160, eos = 280),
      "Frêne commun"     = list(min = 0.15, max = 0.76, sos = 115, peak = 185, eos = 285),
      "Érable sycomore"  = list(min = 0.16, max = 0.79, sos = 100, peak = 170, eos = 290),
      "Peupliers"        = list(min = 0.13, max = 0.81, sos = 95,  peak = 165, eos = 275),
      "Robinier faux-acacia" = list(min = 0.14, max = 0.77, sos = 120, peak = 180, eos = 270),
      "Mélèze d'Europe"  = list(min = 0.10, max = 0.72, sos = 115, peak = 175, eos = 305),
      list(min = 0.15, max = 0.80, sos = 100, peak = 175, eos = 295)  # défaut
    )

    # Ajout variabilité individuelle
    params$sos  <- params$sos + round(runif(1, -10, 10))
    params$peak <- params$peak + round(runif(1, -10, 10))
    params$eos  <- params$eos + round(runif(1, -10, 10))
    params$max  <- params$max + runif(1, -0.05, 0.05)

    # Fonction double logistique
    ndvi <- double_logistic(doy, params$min, params$max,
                            params$sos, params$peak, params$eos)
  }

  # Ajout de bruit
  ndvi <- ndvi + rnorm(length(ndvi), 0, noise_sd)
  ndvi <- pmax(0, pmin(1, ndvi))  # Clamp [0, 1]

  ndvi
}

#' Fonction double logistique pour simuler un cycle phénologique
double_logistic <- function(doy, ndvi_min, ndvi_max, sos, peak, eos) {
  amplitude <- ndvi_max - ndvi_min

  # Montée (verdissement)
  rate_up <- 0.08
  greenup <- 1 / (1 + exp(-rate_up * (doy - sos)))

  # Descente (sénescence)
  rate_down <- 0.06
  senescence <- 1 / (1 + exp(rate_down * (doy - eos)))

  ndvi <- ndvi_min + amplitude * greenup * senescence
  ndvi
}

#' Génération d'un dataset synthétique complet pour les 20 espèces
#' @param n_samples_per_species Nombre d'échantillons par espèce
#' @param year Année de simulation
#' @return data.frame avec colonnes : plot_id, species_code, date, NDVI, + bandes
generate_synthetic_dataset <- function(n_samples_per_species = 50, year = 2021) {
  log_msg("Génération du dataset synthétique TreeSatAI-TS")
  log_msg("  {nrow(SPECIES)} espèces × {n_samples_per_species} échantillons",
          level = "info")

  target_dates <- seq.Date(
    as.Date(paste0(year, "-01-01")),
    as.Date(paste0(year, "-12-31")),
    by = TS_PARAMS$target_interval_days
  )
  n_dates <- length(target_dates)

  all_data <- list()

  for (sp_code in 1:nrow(SPECIES)) {
    sp_name <- SPECIES$french[sp_code]

    for (sample_id in 1:n_samples_per_species) {
      plot_id <- paste0("plot_", sprintf("%02d", sp_code), "_", sprintf("%03d", sample_id))

      set.seed(sp_code * 1000 + sample_id)

      # Simuler NDVI
      ndvi <- simulate_species_ndvi(sp_code, target_dates)

      # Dériver les bandes spectrales à partir du NDVI (simplifié)
      nir <- 0.3 + 0.4 * ndvi + rnorm(n_dates, 0, 0.02)
      red <- nir * (1 - ndvi) / (1 + ndvi + 1e-6) + rnorm(n_dates, 0, 0.01)
      blue <- red * runif(1, 0.7, 0.9) + rnorm(n_dates, 0, 0.01)
      green <- (red + nir) / 3 + rnorm(n_dates, 0, 0.01)
      rededge1 <- (red + nir) / 2 * runif(1, 0.85, 0.95) + rnorm(n_dates, 0, 0.01)
      rededge2 <- (rededge1 + nir) / 2 + rnorm(n_dates, 0, 0.01)
      rededge3 <- nir * runif(1, 0.9, 0.98) + rnorm(n_dates, 0, 0.01)
      nir2 <- nir * runif(1, 0.85, 0.95) + rnorm(n_dates, 0, 0.01)
      swir1 <- 0.2 - 0.1 * ndvi + rnorm(n_dates, 0, 0.02)
      swir2 <- swir1 * runif(1, 0.6, 0.8) + rnorm(n_dates, 0, 0.01)

      # Clamp positif
      bands <- data.frame(
        B02 = pmax(0, blue),
        B03 = pmax(0, green),
        B04 = pmax(0, red),
        B05 = pmax(0, rededge1),
        B06 = pmax(0, rededge2),
        B07 = pmax(0, rededge3),
        B08 = pmax(0, nir),
        B8A = pmax(0, nir2),
        B11 = pmax(0, swir1),
        B12 = pmax(0, swir2)
      )

      # Calculer les indices spectraux à partir des bandes
      idx_ndvi  <- calc_ndvi(bands$B08, bands$B04)
      idx_evi   <- calc_evi(bands$B08, bands$B04, bands$B02)
      idx_ndwi  <- (bands$B03 - bands$B08) / (bands$B03 + bands$B08 + 1e-10)
      idx_nbr   <- (bands$B08 - bands$B12) / (bands$B08 + bands$B12 + 1e-10)
      idx_cri   <- (1 / (bands$B02 + 1e-10)) - (1 / (bands$B03 + 1e-10))
      idx_rendvi <- (bands$B06 - bands$B05) / (bands$B06 + bands$B05 + 1e-10)

      plot_df <- data.frame(
        plot_id      = plot_id,
        species_code = sp_code,
        species_name = sp_name,
        date         = target_dates,
        bands,
        NDVI   = idx_ndvi,
        EVI    = idx_evi,
        NDWI   = idx_ndwi,
        NBR    = idx_nbr,
        CRI    = idx_cri,
        RENDVI = idx_rendvi,
        stringsAsFactors = FALSE
      )

      all_data[[length(all_data) + 1]] <- plot_df
    }
  }

  result <- do.call(rbind, all_data)
  log_msg("Dataset synthétique : {nrow(result)} lignes ({length(unique(result$plot_id))} parcelles)",
          level = "success")

  result
}

#' Conversion du dataset long en features à plat pour la classification
#' @param ts_long Dataset long (sortie de generate_synthetic_dataset)
#' @return data.frame avec 1 ligne par parcelle, toutes features en colonnes
build_feature_matrix <- function(ts_long) {
  log_msg("Construction de la matrice de features...")

  plot_ids <- unique(ts_long$plot_id)
  target_dates <- sort(unique(ts_long$date))
  doy_labels <- format(target_dates, "%j")

  rows <- lapply(plot_ids, function(pid) {
    plot_data <- ts_long[ts_long$plot_id == pid, ]
    plot_data <- plot_data[order(plot_data$date), ]

    species_code <- plot_data$species_code[1]
    species_name <- plot_data$species_name[1]

    # Bandes spectrales comme features temporelles
    band_features <- c()
    bands_ts <- list()

    for (band in S2_BAND_NAMES) {
      if (band %in% names(plot_data)) {
        vals <- plot_data[[band]]
        bands_ts[[band]] <- vals

        # Série temporelle brute
        feat_names <- paste0(band, "_d", doy_labels)
        band_features <- c(band_features, setNames(vals, feat_names))

        # Statistiques temporelles
        stats <- calc_temporal_stats(vals, prefix = band)
        band_features <- c(band_features, stats)
      }
    }

    # Indices spectraux
    idx_features <- c()
    if (length(bands_ts) >= 7) {
      indices <- calc_all_indices(bands_ts)

      for (idx_name in names(indices)) {
        idx_vals <- indices[[idx_name]]
        idx_feat_names <- paste0(idx_name, "_d", doy_labels[seq_along(idx_vals)])
        idx_features <- c(idx_features, setNames(idx_vals, idx_feat_names))

        # Stats
        stats <- calc_temporal_stats(idx_vals, prefix = idx_name)
        idx_features <- c(idx_features, stats)
      }

      # Phénologie sur NDVI lissé
      ndvi_smooth <- smooth_savgol(indices$NDVI)
      pheno <- extract_phenometrics(ndvi_smooth, target_dates)
      idx_features <- c(idx_features, pheno)

      # Fourier
      fourier <- fourier_features(ndvi_smooth, n_harmonics = 3)
      names(fourier) <- paste0("NDVI_", names(fourier))
      idx_features <- c(idx_features, fourier)
    }

    c(plot_id = pid, species_code = species_code, species_name = species_name,
      band_features, idx_features)
  })

  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)

  # Conversion numérique
  numeric_cols <- setdiff(names(df), c("plot_id", "species_name"))
  df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

  log_msg("Matrice de features : {nrow(df)} parcelles × {ncol(df)} colonnes",
          level = "success")

  df
}

