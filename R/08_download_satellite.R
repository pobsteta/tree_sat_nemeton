#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Téléchargement Sentinel-2 & Sentinel-1
# via rstac + Microsoft Planetary Computer
#
# Sentinel-2 L2A : 10 bandes optiques, résolution 10-20m
# Sentinel-1 RTC : rétrodiffusion radar VV/VH, déjà calibré, résolution 10m
#
# AUCUNE AUTHENTIFICATION requise (Planetary Computer est gratuit et ouvert)
#
# Les données sont des Cloud-Optimized GeoTIFF (COG) : seule l'emprise de
# l'AOI est téléchargée via /vsicurl/, ce qui est très efficace.
# ==============================================================================

source(file.path(here::here(), "R", "00_config.R"))
source(file.path(here::here(), "R", "01_utils.R"))

# ==============================================================================
# CONFIGURATION — MICROSOFT PLANETARY COMPUTER
# ==============================================================================

STAC_CONFIG <- list(
  stac_url      = "https://planetarycomputer.microsoft.com/api/stac/v1",
  s2_collection = "sentinel-2-l2a",
  s1_collection = "sentinel-1-rtc",

  # Bandes S2 à télécharger
  s2_bands = c("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"),

  # Assets S1 (Planetary Computer utilise des noms en minuscules)
  s1_assets = c("vv", "vh"),

  # Limites
  max_results = 500,
  retry_wait  = c(2, 4, 8, 16)
)

# Bandes Sentinel-1
S1_BANDS <- list(
  VV = list(name = "VV", polarization = "VV", description = "Co-polarisation verticale"),
  VH = list(name = "VH", polarization = "VH", description = "Cross-polarisation")
)

S1_BAND_NAMES <- names(S1_BANDS)

# ==============================================================================
# 1. RECHERCHE DE PRODUITS VIA rstac
# ==============================================================================

#' Recherche de produits Sentinel-2 L2A via Planetary Computer
#' @param aoi sf object — zone d'intérêt
#' @param start_date Date de début "YYYY-MM-DD"
#' @param end_date Date de fin "YYYY-MM-DD"
#' @param max_cloud Couverture nuageuse max (%)
#' @return Liste avec data.frame scenes + objet rstac items signé
search_sentinel2 <- function(aoi, start_date, end_date, max_cloud = 30) {
  log_msg("Recherche Sentinel-2 L2A via Planetary Computer...")

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  bbox <- as.numeric(sf::st_bbox(aoi_wgs84))

  log_msg("  Bbox : [{round(bbox[1],4)}, {round(bbox[2],4)}, {round(bbox[3],4)}, {round(bbox[4],4)}]")
  log_msg("  Période : {start_date} → {end_date}, nuages ≤ {max_cloud}%")

  datetime_str <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")

  # Planetary Computer ne supporte pas ext_query ; on utilise get_request()
  # puis on filtre côté client par eo:cloud_cover
  items <- tryCatch({
    rstac::stac(STAC_CONFIG$stac_url) |>
      rstac::stac_search(
        collections = STAC_CONFIG$s2_collection,
        bbox        = bbox,
        datetime    = datetime_str,
        limit       = STAC_CONFIG$max_results
      ) |>
      rstac::get_request()
  }, error = function(e) {
    cli::cli_alert_danger("Erreur STAC S2 : {e$message}")
    return(NULL)
  })

  if (is.null(items) || length(items$features) == 0) {
    log_msg("  Aucune scène S2 trouvée", level = "warning")
    return(NULL)
  }

  # Filtrage côté client : garder uniquement les scènes ≤ max_cloud
  items$features <- Filter(function(feat) {
    cc <- feat$properties$`eo:cloud_cover`
    !is.null(cc) && cc <= max_cloud
  }, items$features)

  if (length(items$features) == 0) {
    log_msg("  Aucune scène S2 ≤ {max_cloud}% de nuages", level = "warning")
    return(NULL)
  }

  # Signer les URLs (Planetary Computer — gratuit, pas de compte)
  items_signed <- rstac::items_sign(items, sign_fn = rstac::sign_planetary_computer())

  # Extraire les métadonnées
  scenes <- lapply(items_signed$features, function(feat) {
    data.frame(
      id           = feat$id %||% NA_character_,
      datetime     = feat$properties$datetime %||% NA_character_,
      date         = as.Date(substr(feat$properties$datetime %||% "", 1, 10)),
      cloud_cover  = feat$properties$`eo:cloud_cover` %||% NA_real_,
      platform     = feat$properties$platform %||% NA_character_,
      product_type = "S2_L2A",
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, scenes)
  result <- result[order(result$date), ]

  log_msg("  {nrow(result)} scènes S2 L2A trouvées ({min(result$date)} → {max(result$date)})",
          level = "success")

  result$month <- format(result$date, "%Y-%m")
  monthly <- table(result$month)
  log_msg("  Répartition mensuelle : {paste(names(monthly), monthly, sep=':', collapse=', ')}")

  # Retourner data.frame + items signés
  list(scenes = result, items = items_signed)
}

#' Recherche de produits Sentinel-1 RTC via Planetary Computer
#' @param aoi sf object — zone d'intérêt
#' @param start_date Date de début
#' @param end_date Date de fin
#' @param orbit_direction Direction d'orbite ("ascending", "descending", ou NULL)
#' @return Liste avec data.frame scenes + objet rstac items signé
search_sentinel1 <- function(aoi, start_date, end_date, orbit_direction = NULL) {
  log_msg("Recherche Sentinel-1 RTC via Planetary Computer...")

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  bbox <- as.numeric(sf::st_bbox(aoi_wgs84))

  log_msg("  Bbox : [{round(bbox[1],4)}, {round(bbox[2],4)}, {round(bbox[3],4)}, {round(bbox[4],4)}]")
  log_msg("  Période : {start_date} → {end_date}")

  datetime_str <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")

  # Planetary Computer : get_request() + filtrage côté client
  items <- tryCatch({
    rstac::stac(STAC_CONFIG$stac_url) |>
      rstac::stac_search(
        collections = STAC_CONFIG$s1_collection,
        bbox        = bbox,
        datetime    = datetime_str,
        limit       = STAC_CONFIG$max_results
      ) |>
      rstac::get_request()
  }, error = function(e) {
    cli::cli_alert_danger("Erreur STAC S1 : {e$message}")
    return(NULL)
  })

  if (is.null(items) || length(items$features) == 0) {
    log_msg("  Aucun produit S1 trouvé", level = "warning")
    return(NULL)
  }

  # Filtrage côté client : direction d'orbite si spécifiée
  if (!is.null(orbit_direction)) {
    target <- tolower(orbit_direction)
    items$features <- Filter(function(feat) {
      orb <- feat$properties$`sat:orbit_state`
      !is.null(orb) && tolower(orb) == target
    }, items$features)

    if (length(items$features) == 0) {
      log_msg("  Aucun produit S1 en orbite {orbit_direction}", level = "warning")
      return(NULL)
    }
  }

  items_signed <- rstac::items_sign(items, sign_fn = rstac::sign_planetary_computer())

  scenes <- lapply(items_signed$features, function(feat) {
    data.frame(
      id           = feat$id %||% NA_character_,
      name         = feat$id %||% NA_character_,
      datetime     = feat$properties$datetime %||% NA_character_,
      date         = as.Date(substr(feat$properties$datetime %||% "", 1, 10)),
      platform     = feat$properties$platform %||% NA_character_,
      product_type = "S1_RTC",
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, scenes)
  result <- result[order(result$date), ]

  log_msg("  {nrow(result)} produits S1 RTC trouvés", level = "success")

  list(scenes = result, items = items_signed)
}

# ==============================================================================
# 2. TÉLÉCHARGEMENT VIA /vsicurl/ (COG streaming)
# ==============================================================================

#' Lecture d'un raster distant (COG) et crop sur l'AOI
#' Utilise /vsicurl/ pour ne télécharger que l'emprise nécessaire.
#' Fallback sur download.file si vsicurl échoue.
#'
#' @param url URL du COG distant (signée)
#' @param aoi_vect SpatVector de l'AOI (dans le CRS du raster ou WGS84)
#' @param buffer_m Marge autour de l'AOI (mètres)
#' @return SpatRaster croppé, ou NULL en cas d'erreur
read_remote_cog <- function(url, aoi_vect, buffer_m = 1000) {
  # Tentative 1 : /vsicurl/ (efficace, lecture partielle du COG)
  r <- tryCatch({
    r <- terra::rast(paste0("/vsicurl/", url))
    aoi_native <- terra::project(aoi_vect, terra::crs(r))
    crop_ext <- terra::ext(terra::buffer(aoi_native, buffer_m))
    terra::crop(r, crop_ext)
  }, error = function(e) NULL)

  if (!is.null(r)) return(r)

  # Tentative 2 : téléchargement complet + lecture locale
  r <- tryCatch({
    tmp <- tempfile(fileext = ".tif")
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    r <- terra::rast(tmp)
    aoi_native <- terra::project(aoi_vect, terra::crs(r))
    crop_ext <- terra::ext(terra::buffer(aoi_native, buffer_m))
    r <- terra::crop(r, crop_ext)
    unlink(tmp)
    r
  }, error = function(e) NULL)

  r
}

#' Téléchargement des bandes S2 pour les scènes sélectionnées
#' @param items_signed Objet rstac items signé
#' @param selected_ids Vecteur d'IDs de scènes sélectionnées
#' @param aoi sf object
#' @param bands_dir Répertoire de sortie pour les bandes
#' @return Nombre de bandes téléchargées
download_s2_bands <- function(items_signed, selected_ids, aoi, bands_dir) {
  dir.create(bands_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  aoi_vect  <- terra::vect(aoi_wgs84)

  # Filtrer les items sélectionnés
  sel_features <- items_signed$features[
    sapply(items_signed$features, function(f) f$id %in% selected_ids)
  ]

  n_total <- length(sel_features) * length(STAC_CONFIG$s2_bands)
  pb <- cli::cli_progress_bar("Téléchargement bandes S2", total = n_total)
  n_ok <- 0

  for (feat in sel_features) {
    date_str <- gsub("-", "", substr(feat$properties$datetime, 1, 10))

    for (band in STAC_CONFIG$s2_bands) {
      out_path <- file.path(bands_dir, paste0(date_str, "_", band, ".tif"))

      if (file.exists(out_path)) {
        n_ok <- n_ok + 1
        cli::cli_progress_update(id = pb)
        next
      }

      url <- feat$assets[[band]]$href
      if (is.null(url)) {
        cli::cli_progress_update(id = pb)
        next
      }

      r <- read_remote_cog(url, aoi_vect)
      if (!is.null(r)) {
        terra::writeRaster(r, out_path, overwrite = TRUE)
        n_ok <- n_ok + 1
      } else {
        cli::cli_alert_warning("  Échec : {date_str}_{band}")
      }

      cli::cli_progress_update(id = pb)
    }
  }

  cli::cli_progress_done(id = pb)
  log_msg("  {n_ok}/{n_total} bandes S2 téléchargées dans {bands_dir}", level = "success")

  n_ok
}

#' Téléchargement des bandes S1 (VV, VH) avec conversion en dB
#' @param items_signed Objet rstac items signé
#' @param selected_ids Vecteur d'IDs de scènes sélectionnées
#' @param aoi sf object
#' @param bands_dir Répertoire de sortie
#' @return Nombre de bandes téléchargées
download_s1_bands <- function(items_signed, selected_ids, aoi, bands_dir) {
  dir.create(bands_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  aoi_vect  <- terra::vect(aoi_wgs84)

  sel_features <- items_signed$features[
    sapply(items_signed$features, function(f) f$id %in% selected_ids)
  ]

  n_total <- length(sel_features) * length(STAC_CONFIG$s1_assets)
  pb <- cli::cli_progress_bar("Téléchargement bandes S1", total = n_total)
  n_ok <- 0

  for (feat in sel_features) {
    date_str <- gsub("-", "", substr(feat$properties$datetime, 1, 10))

    for (asset_name in STAC_CONFIG$s1_assets) {
      pol <- toupper(asset_name)  # vv → VV
      out_path <- file.path(bands_dir, paste0(date_str, "_", pol, ".tif"))

      if (file.exists(out_path)) {
        n_ok <- n_ok + 1
        cli::cli_progress_update(id = pb)
        next
      }

      url <- feat$assets[[asset_name]]$href
      if (is.null(url)) {
        cli::cli_progress_update(id = pb)
        next
      }

      r <- read_remote_cog(url, aoi_vect)
      if (!is.null(r)) {
        # Conversion gamma0 linéaire → dB
        vals <- terra::values(r)
        vals[vals <= 0] <- NA
        vals <- 10 * log10(vals)
        terra::values(r) <- vals

        terra::writeRaster(r, out_path, overwrite = TRUE)
        n_ok <- n_ok + 1
      } else {
        cli::cli_alert_warning("  Échec : {date_str}_{pol}")
      }

      cli::cli_progress_update(id = pb)
    }
  }

  cli::cli_progress_done(id = pb)
  log_msg("  {n_ok}/{n_total} bandes S1 téléchargées dans {bands_dir}", level = "success")

  n_ok
}

# ==============================================================================
# 3. PIPELINES DE TÉLÉCHARGEMENT
# ==============================================================================

#' Téléchargement Sentinel-2 pour une AOI et une année
#' @param aoi sf object
#' @param year Année
#' @param max_cloud Couverture nuageuse max
#' @param output_dir Répertoire de sortie
#' @param max_scenes Nombre max de scènes
#' @return Chemin vers le répertoire des données
download_s2_for_aoi <- function(aoi, year = 2021, max_cloud = 30,
                                 output_dir = file.path(RAW_DIR, "sentinel2"),
                                 max_scenes = NULL) {
  cli::cli_h2("Téléchargement Sentinel-2 L2A (Planetary Computer)")

  start_date <- paste0(year, "-01-01")
  end_date   <- paste0(year, "-12-31")

  search_result <- search_sentinel2(aoi, start_date, end_date, max_cloud)
  if (is.null(search_result)) return(NULL)

  scenes      <- search_result$scenes
  items_signed <- search_result$items

  # Sélection temporelle : 1 scène tous les ~10 jours
  scenes <- select_best_scenes(scenes, interval_days = 10)
  log_msg("  {nrow(scenes)} scènes sélectionnées (1 / 10 jours)")

  if (!is.null(max_scenes)) {
    scenes <- scenes[1:min(max_scenes, nrow(scenes)), ]
    log_msg("  Limité à {nrow(scenes)} scènes (max_scenes)")
  }

  if (interactive()) {
    confirm <- readline(glue::glue(
      "Télécharger {nrow(scenes)} scènes S2 × {length(STAC_CONFIG$s2_bands)} bandes ? (o/n) : "
    ))
    if (!tolower(confirm) %in% c("o", "oui", "y", "yes")) {
      log_msg("Téléchargement annulé", level = "warning")
      return(NULL)
    }
  }

  s2_dir    <- file.path(output_dir, paste0("S2_L2A_", year))
  bands_dir <- file.path(s2_dir, "bands")

  download_s2_bands(items_signed, scenes$id, aoi, bands_dir)

  s2_dir
}

#' Téléchargement Sentinel-1 pour une AOI et une année
#' @param aoi sf object
#' @param year Année
#' @param output_dir Répertoire de sortie
#' @param max_scenes Nombre max de scènes
#' @return Chemin vers le répertoire des données
download_s1_for_aoi <- function(aoi, year = 2021,
                                 output_dir = file.path(RAW_DIR, "sentinel1"),
                                 max_scenes = NULL) {
  cli::cli_h2("Téléchargement Sentinel-1 RTC (Planetary Computer)")

  start_date <- paste0(year, "-01-01")
  end_date   <- paste0(year, "-12-31")

  search_result <- search_sentinel1(aoi, start_date, end_date, orbit_direction = "descending")
  if (is.null(search_result)) return(NULL)

  scenes       <- search_result$scenes
  items_signed <- search_result$items

  scenes <- select_best_scenes_s1(scenes, interval_days = 12)
  log_msg("  {nrow(scenes)} scènes S1 sélectionnées")

  if (!is.null(max_scenes)) {
    scenes <- scenes[1:min(max_scenes, nrow(scenes)), ]
  }

  if (interactive()) {
    confirm <- readline(glue::glue(
      "Télécharger {nrow(scenes)} scènes S1 (VV + VH) ? (o/n) : "
    ))
    if (!tolower(confirm) %in% c("o", "oui", "y", "yes")) {
      log_msg("Téléchargement S1 annulé", level = "warning")
      return(NULL)
    }
  }

  s1_dir    <- file.path(output_dir, paste0("S1_RTC_", year))
  bands_dir <- file.path(s1_dir, "bands")

  download_s1_bands(items_signed, scenes$id, aoi, bands_dir)

  s1_dir
}

# ==============================================================================
# 4. PIPELINE COMPLET S2 + S1
# ==============================================================================

#' Téléchargement automatique de toutes les données satellite pour une AOI
#' Aucune authentification requise (Planetary Computer est gratuit)
#'
#' @param aoi_path Chemin vers l'AOI (GeoPackage, Shapefile...)
#' @param year Année d'analyse
#' @param download_s2 Télécharger Sentinel-2 (défaut TRUE)
#' @param download_s1 Télécharger Sentinel-1 (défaut TRUE)
#' @param max_cloud Couverture nuageuse max S2 (%)
#' @param output_dir Répertoire de sortie
#' @return Liste avec les chemins vers les données téléchargées
download_satellite_data <- function(aoi_path, year = 2021,
                                     download_s2 = TRUE, download_s1 = TRUE,
                                     max_cloud = 30,
                                     output_dir = RAW_DIR) {
  cli::cli_h1("Téléchargement des données satellite (Planetary Computer)")
  cli::cli_text("Aucune authentification requise")
  t_start <- Sys.time()

  aoi <- sf::st_read(aoi_path, quiet = TRUE)
  if (nrow(aoi) > 1) aoi <- sf::st_union(aoi) |> sf::st_as_sf()

  aoi_area_ha <- as.numeric(sf::st_area(sf::st_transform(aoi, 2154))) / 10000
  log_msg("AOI : {round(aoi_area_ha, 1)} ha")

  results <- list(s2_dir = NULL, s1_dir = NULL)

  if (download_s2) {
    results$s2_dir <- download_s2_for_aoi(
      aoi, year = year, max_cloud = max_cloud,
      output_dir = output_dir
    )
  }

  if (download_s1) {
    results$s1_dir <- download_s1_for_aoi(
      aoi, year = year,
      output_dir = output_dir
    )
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  cli::cli_h2("Téléchargement terminé ({round(t_elapsed, 1)} min)")

  if (!is.null(results$s2_dir)) {
    n_s2 <- length(list.files(file.path(results$s2_dir, "bands"),
                              pattern = "\\.tif$", recursive = TRUE))
    cli::cli_alert_success("Sentinel-2 : {n_s2} bandes → {results$s2_dir}")
  }
  if (!is.null(results$s1_dir)) {
    n_s1 <- length(list.files(file.path(results$s1_dir, "bands"),
                              pattern = "\\.tif$", recursive = TRUE))
    cli::cli_alert_success("Sentinel-1 : {n_s1} bandes → {results$s1_dir}")
  }

  invisible(results)
}

# ==============================================================================
# 5. FONCTIONS UTILITAIRES
# ==============================================================================

#' Sélection des meilleures scènes S2 (1 par intervalle, moins de nuages)
select_best_scenes <- function(scenes, interval_days = 10) {
  if (nrow(scenes) == 0) return(scenes)

  scenes <- scenes[order(scenes$date, scenes$cloud_cover), ]

  selected <- scenes[1, ]
  last_date <- scenes$date[1]

  for (i in 2:nrow(scenes)) {
    days_since <- as.numeric(difftime(scenes$date[i], last_date, units = "days"))
    if (days_since >= interval_days) {
      selected <- rbind(selected, scenes[i, ])
      last_date <- scenes$date[i]
    }
  }

  selected
}

#' Sélection des meilleures scènes S1
select_best_scenes_s1 <- function(scenes, interval_days = 12) {
  if (nrow(scenes) == 0) return(scenes)

  scenes <- scenes[order(scenes$date), ]

  selected <- scenes[1, ]
  last_date <- scenes$date[1]

  for (i in 2:nrow(scenes)) {
    days_since <- as.numeric(difftime(scenes$date[i], last_date, units = "days"))
    if (days_since >= interval_days) {
      selected <- rbind(selected, scenes[i, ])
      last_date <- scenes$date[i]
    }
  }

  selected
}

#' Construction du cube Sentinel-1 (VV + VH time series)
#' @param s1_bands_dir Répertoire des bandes S1 traitées
#' @param aoi sf object
#' @param year Année
#' @param resolution Résolution cible
#' @return Liste de SpatRaster
build_s1_cube <- function(s1_bands_dir, aoi, year = 2021, resolution = 10) {
  log_msg("Construction du cube Sentinel-1")

  aoi_proj <- sf::st_transform(aoi, 2154)
  aoi_vect <- terra::vect(aoi_proj)
  aoi_ext  <- terra::ext(aoi_vect)

  files <- list.files(s1_bands_dir, pattern = "\\.tif$", full.names = TRUE)
  if (length(files) == 0) {
    log_msg("  Aucun fichier S1 trouvé", level = "warning")
    return(NULL)
  }

  file_info <- data.frame(
    path = files,
    date = as.Date(stringr::str_extract(basename(files), "\\d{8}"), format = "%Y%m%d"),
    pol  = stringr::str_extract(basename(files), "VV|VH"),
    stringsAsFactors = FALSE
  )

  file_info <- file_info[format(file_info$date, "%Y") == as.character(year), ]
  dates <- sort(unique(file_info$date))

  log_msg("  {length(dates)} dates S1 en {year}")

  cube_list <- list()

  for (d in dates) {
    date_files <- file_info[file_info$date == d, ]
    stack <- list()

    for (pol in c("VV", "VH")) {
      pf <- date_files[date_files$pol == pol, ]
      if (nrow(pf) > 0) {
        r <- terra::rast(pf$path[1])
        if (!terra::same.crs(r, terra::crs("EPSG:2154"))) {
          r <- terra::project(r, "EPSG:2154", method = "bilinear")
        }
        r <- terra::crop(r, aoi_ext)
        r <- terra::mask(r, aoi_vect)
        names(r) <- paste0(pol, "_", format(d, "%Y%m%d"))
        stack[[pol]] <- r
      }
    }

    # Ratio VV/VH (en dB = soustraction)
    if ("VV" %in% names(stack) && "VH" %in% names(stack)) {
      ratio <- stack$VV - stack$VH
      names(ratio) <- paste0("VV_VH_ratio_", format(d, "%Y%m%d"))
      stack$ratio <- ratio
    }

    if (length(stack) >= 2) {
      cube_list[[as.character(d)]] <- terra::rast(stack)
    }
  }

  log_msg("  Cube S1 : {length(cube_list)} dates × {length(S1_BAND_NAMES)} polarisations + ratio",
          level = "success")

  cube_list
}

# ==============================================================================
# 6. EXTRACTION DE FEATURES S1
# ==============================================================================

#' Calcul d'indices radar pour la classification forestière
calc_radar_indices <- function(vv, vh) {
  list(
    VV_VH_ratio = vv - vh,
    RVI = 4 * 10^(vh/10) / (10^(vv/10) + 10^(vh/10)),
    RFDI = (10^(vv/10) - 10^(vh/10)) /
           (10^(vv/10) + 10^(vh/10))
  )
}

#' Statistiques temporelles Sentinel-1 pour un pixel
calc_s1_temporal_features <- function(vv_ts, vh_ts) {
  features <- c()

  features <- c(features, calc_temporal_stats(vv_ts, prefix = "S1_VV"))
  features <- c(features, calc_temporal_stats(vh_ts, prefix = "S1_VH"))

  ratio_ts <- vv_ts - vh_ts
  features <- c(features, calc_temporal_stats(ratio_ts, prefix = "S1_ratio"))

  rvi_ts <- 4 * 10^(vh_ts/10) / (10^(vv_ts/10) + 10^(vh_ts/10))
  features <- c(features, calc_temporal_stats(rvi_ts, prefix = "S1_RVI"))

  n <- length(vv_ts)
  if (n >= 4) {
    q1 <- 1:floor(n/4)
    q2 <- (floor(n/4)+1):(floor(n/2))
    q3 <- (floor(n/2)+1):(floor(3*n/4))
    q4 <- (floor(3*n/4)+1):n

    features <- c(features,
      S1_VH_q1 = mean(vh_ts[q1], na.rm = TRUE),
      S1_VH_q2 = mean(vh_ts[q2], na.rm = TRUE),
      S1_VH_q3 = mean(vh_ts[q3], na.rm = TRUE),
      S1_VH_q4 = mean(vh_ts[q4], na.rm = TRUE),
      S1_VH_summer_winter = mean(vh_ts[q2], na.rm = TRUE) - mean(vh_ts[q4], na.rm = TRUE)
    )
  }

  features
}

# ==============================================================================
# 7. NULL-SAFE OPERATOR
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

cli::cli_alert_success("Module satellite chargé (Planetary Computer — sans authentification)")
cli::cli_text("Utilisation :")
cli::cli_text('  {.code download_satellite_data("aoi.gpkg", year = 2023)}')
cli::cli_text('  {.code download_s2_for_aoi(aoi, year = 2023)}')
cli::cli_text('  {.code download_s1_for_aoi(aoi, year = 2023)}')
