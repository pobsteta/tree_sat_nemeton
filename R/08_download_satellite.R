#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Téléchargement Sentinel-2 & Sentinel-1
# via rstac + Copernicus Data Space Ecosystem (CDSE)
#
# Sentinel-2 L2A : 10 bandes optiques, résolution 10-20m
# Sentinel-1 GRD : rétrodiffusion radar VV/VH, résolution 10m
#
# Prérequis : compte gratuit sur https://dataspace.copernicus.eu
# ==============================================================================

source(file.path(here::here(), "R", "00_config.R"))
source(file.path(here::here(), "R", "01_utils.R"))

# ==============================================================================
# CONFIGURATION CDSE
# ==============================================================================

CDSE_CONFIG <- list(
  stac_url     = "https://catalogue.dataspace.copernicus.eu/stac",
  token_url    = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
  download_url = "https://zipper.dataspace.copernicus.eu/odata/v1",

  # Collections STAC
  s2_collection = "sentinel-2-l2a",
  s1_collection = "SENTINEL-1",

  # Limites
  max_results = 1000,
  retry_max   = 4,
  retry_wait  = c(2, 4, 8, 16)
)

# Bandes Sentinel-1
S1_BANDS <- list(
  VV = list(name = "VV", polarization = "VV", description = "Co-polarisation verticale"),
  VH = list(name = "VH", polarization = "VH", description = "Cross-polarisation")
)

S1_BAND_NAMES <- names(S1_BANDS)

# ==============================================================================
# 1. AUTHENTIFICATION CDSE
# ==============================================================================

#' Obtention d'un token d'accès CDSE via OAuth2
#' @param username Email du compte CDSE
#' @param password Mot de passe CDSE
#' @return Token d'accès (chaîne de caractères)
cdse_get_token <- function(username = NULL, password = NULL) {
  if (is.null(username)) username <- Sys.getenv("CDSE_USERNAME", unset = NA)
  if (is.null(password)) password <- Sys.getenv("CDSE_PASSWORD", unset = NA)

  if (is.na(username) || is.na(password)) {
    if (interactive()) {
      cli::cli_h3("Authentification Copernicus Data Space")
      cli::cli_text("Compte gratuit : https://dataspace.copernicus.eu")
      cli::cli_text("(Ou définir CDSE_USERNAME / CDSE_PASSWORD dans .Renviron)")
      cli::cli_text("")
      username <- readline("Email CDSE : ")
      password <- readline("Mot de passe : ")
    } else {
      cli::cli_alert_danger("Identifiants CDSE manquants.")
      cli::cli_text("Définir dans ~/.Renviron :")
      cli::cli_text("  CDSE_USERNAME=votre@email.com")
      cli::cli_text("  CDSE_PASSWORD=motdepasse")
      return(NULL)
    }
  }

  log_msg("Authentification CDSE...")

  resp <- tryCatch({
    httr2::request(CDSE_CONFIG$token_url) |>
      httr2::req_body_form(
        grant_type = "password",
        username   = username,
        password   = password,
        client_id  = "cdse-public"
      ) |>
      httr2::req_timeout(15) |>
      httr2::req_perform()
  }, error = function(e) {
    cli::cli_alert_danger("Erreur d'authentification : {e$message}")
    return(NULL)
  })

  if (is.null(resp)) return(NULL)

  token_data <- httr2::resp_body_json(resp)
  token <- token_data$access_token

  if (is.null(token)) {
    cli::cli_alert_danger("Token non obtenu — vérifiez vos identifiants")
    return(NULL)
  }

  log_msg("  Authentifié (token valide {token_data$expires_in}s)", level = "success")
  token
}

# ==============================================================================
# 2. RECHERCHE DE PRODUITS VIA rstac
# ==============================================================================

#' Recherche de produits Sentinel-2 L2A sur une AOI via rstac
#' @param aoi sf object — zone d'intérêt
#' @param start_date Date de début "YYYY-MM-DD"
#' @param end_date Date de fin "YYYY-MM-DD"
#' @param max_cloud Couverture nuageuse max (%)
#' @return data.frame avec métadonnées des produits
search_sentinel2 <- function(aoi, start_date, end_date, max_cloud = 30) {
  log_msg("Recherche Sentinel-2 L2A via rstac...")

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  bbox <- as.numeric(sf::st_bbox(aoi_wgs84))

  log_msg("  Bbox : [{round(bbox[1],4)}, {round(bbox[2],4)}, {round(bbox[3],4)}, {round(bbox[4],4)}]")
  log_msg("  Période : {start_date} → {end_date}, nuages ≤ {max_cloud}%")

  datetime_str <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")

  # Requête STAC via rstac
  items <- tryCatch({
    rstac::stac(CDSE_CONFIG$stac_url) |>
      rstac::stac_search(
        collections = CDSE_CONFIG$s2_collection,
        bbox        = bbox,
        datetime    = datetime_str,
        limit       = CDSE_CONFIG$max_results
      ) |>
      rstac::ext_query(`eo:cloud_cover` = list(`lte` = max_cloud)) |>
      rstac::post_request()
  }, error = function(e) {
    cli::cli_alert_danger("Erreur STAC S2 : {e$message}")
    return(NULL)
  })

  if (is.null(items) || length(items$features) == 0) {
    log_msg("  Aucune scène S2 trouvée", level = "warning")
    return(NULL)
  }

  # Pagination : récupérer toutes les pages
  all_items <- items
  while (rstac::items_length(all_items) < rstac::items_matched(items) &&
         !is.null(tryCatch(rstac::items_next(all_items), error = function(e) NULL))) {
    next_page <- tryCatch(
      rstac::items_next(all_items) |> rstac::post_request(),
      error = function(e) NULL
    )
    if (is.null(next_page) || length(next_page$features) == 0) break
    all_items$features <- c(all_items$features, next_page$features)
  }

  # Extraire les métadonnées
  scenes <- lapply(all_items$features, function(feat) {
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

  result
}

#' Recherche de produits Sentinel-1 GRD sur une AOI via rstac
#' @param aoi sf object — zone d'intérêt
#' @param start_date Date de début
#' @param end_date Date de fin
#' @param orbit_direction Direction d'orbite ("ASCENDING", "DESCENDING", ou NULL)
#' @return data.frame avec métadonnées des produits
search_sentinel1 <- function(aoi, start_date, end_date, orbit_direction = NULL) {
  log_msg("Recherche Sentinel-1 GRD via rstac...")

  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  bbox <- as.numeric(sf::st_bbox(aoi_wgs84))

  log_msg("  Bbox : [{round(bbox[1],4)}, {round(bbox[2],4)}, {round(bbox[3],4)}, {round(bbox[4],4)}]")
  log_msg("  Période : {start_date} → {end_date}")

  datetime_str <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")

  # Requête STAC via rstac
  q <- rstac::stac(CDSE_CONFIG$stac_url) |>
    rstac::stac_search(
      collections = CDSE_CONFIG$s1_collection,
      bbox        = bbox,
      datetime    = datetime_str,
      limit       = CDSE_CONFIG$max_results
    )

  # Filtre direction d'orbite si spécifiée
  if (!is.null(orbit_direction)) {
    q <- q |> rstac::ext_query(
      `sat:orbit_state` = list(`eq` = tolower(orbit_direction))
    )
  }

  items <- tryCatch({
    q |> rstac::post_request()
  }, error = function(e) {
    cli::cli_alert_danger("Erreur STAC S1 : {e$message}")
    return(NULL)
  })

  if (is.null(items) || length(items$features) == 0) {
    log_msg("  Aucun produit S1 trouvé", level = "warning")
    return(NULL)
  }

  scenes <- lapply(items$features, function(feat) {
    data.frame(
      id           = feat$id %||% NA_character_,
      name         = feat$id %||% NA_character_,
      datetime     = feat$properties$datetime %||% NA_character_,
      date         = as.Date(substr(feat$properties$datetime %||% "", 1, 10)),
      platform     = feat$properties$platform %||% NA_character_,
      product_type = "S1_GRD",
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, scenes)
  result <- result[order(result$date), ]

  log_msg("  {nrow(result)} produits S1 GRD trouvés", level = "success")

  result
}

# ==============================================================================
# 3. TÉLÉCHARGEMENT
# ==============================================================================

#' Téléchargement d'un produit Sentinel depuis CDSE
#' @param product_id ID du produit
#' @param product_name Nom du produit (pour le nom de fichier)
#' @param token Token d'accès CDSE
#' @param output_dir Répertoire de destination
#' @return Chemin vers le fichier téléchargé
download_product <- function(product_id, product_name, token, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  dest_file <- file.path(output_dir, paste0(product_name, ".zip"))

  if (file.exists(dest_file)) {
    log_msg("  Déjà téléchargé : {basename(dest_file)}")
    return(dest_file)
  }

  download_url <- glue::glue("{CDSE_CONFIG$download_url}/Products({product_id})/$value")

  for (attempt in seq_along(CDSE_CONFIG$retry_wait)) {
    result <- tryCatch({
      httr2::request(download_url) |>
        httr2::req_auth_bearer_token(token) |>
        httr2::req_timeout(600) |>
        httr2::req_perform(path = dest_file)
      TRUE
    }, error = function(e) {
      cli::cli_alert_warning("  Tentative {attempt} échouée : {e$message}")
      FALSE
    })

    if (result) {
      size_mb <- round(file.info(dest_file)$size / 1e6, 1)
      log_msg("  Téléchargé : {basename(dest_file)} ({size_mb} Mo)", level = "success")
      return(dest_file)
    }

    if (attempt < length(CDSE_CONFIG$retry_wait)) {
      wait_s <- CDSE_CONFIG$retry_wait[attempt]
      log_msg("  Attente {wait_s}s avant nouvelle tentative...")
      Sys.sleep(wait_s)

      token <- cdse_get_token()
      if (is.null(token)) return(NULL)
    }
  }

  cli::cli_alert_danger("Échec du téléchargement après {length(CDSE_CONFIG$retry_wait)} tentatives")
  return(NULL)
}

#' Téléchargement de tous les produits S2 pour une AOI et une année
#' @param aoi sf object
#' @param year Année
#' @param max_cloud Couverture nuageuse max
#' @param output_dir Répertoire de sortie
#' @param token Token CDSE (NULL = demande interactive)
#' @param max_scenes Nombre max de scènes
#' @return Chemin vers le répertoire des données téléchargées
download_s2_for_aoi <- function(aoi, year = 2021, max_cloud = 30,
                                 output_dir = file.path(RAW_DIR, "sentinel2"),
                                 token = NULL, max_scenes = NULL) {
  cli::cli_h2("Téléchargement Sentinel-2 L2A")

  if (is.null(token)) token <- cdse_get_token()
  if (is.null(token)) return(NULL)

  start_date <- paste0(year, "-01-01")
  end_date   <- paste0(year, "-12-31")
  scenes <- search_sentinel2(aoi, start_date, end_date, max_cloud)
  if (is.null(scenes)) return(NULL)

  # Sélection temporelle : 1 scène tous les ~10 jours
  scenes <- select_best_scenes(scenes, interval_days = 10)
  log_msg("  {nrow(scenes)} scènes sélectionnées (1 / 10 jours)")

  if (!is.null(max_scenes)) {
    scenes <- scenes[1:min(max_scenes, nrow(scenes)), ]
    log_msg("  Limité à {nrow(scenes)} scènes (max_scenes)")
  }

  est_size_gb <- nrow(scenes) * 0.8
  log_msg("  Taille estimée : ~{round(est_size_gb, 1)} Go")

  if (interactive()) {
    confirm <- readline(glue::glue(
      "Télécharger {nrow(scenes)} scènes S2 (~{round(est_size_gb, 1)} Go) ? (o/n) : "
    ))
    if (!tolower(confirm) %in% c("o", "oui", "y", "yes")) {
      log_msg("Téléchargement annulé", level = "warning")
      return(NULL)
    }
  }

  s2_dir <- file.path(output_dir, paste0("S2_L2A_", year))
  dir.create(s2_dir, showWarnings = FALSE, recursive = TRUE)

  pb <- cli::cli_progress_bar("Téléchargement S2", total = nrow(scenes))

  downloaded <- c()
  for (i in seq_len(nrow(scenes))) {
    scene <- scenes[i, ]
    log_msg("  [{i}/{nrow(scenes)}] {scene$id} ({scene$date}, {round(scene$cloud_cover)}% nuages)")

    zip_path <- download_product(scene$id, scene$id, token, s2_dir)

    if (!is.null(zip_path)) {
      extract_dir <- file.path(s2_dir, tools::file_path_sans_ext(basename(zip_path)))
      if (!dir.exists(extract_dir)) {
        unzip(zip_path, exdir = s2_dir)
      }
      downloaded <- c(downloaded, zip_path)
    }

    cli::cli_progress_update(id = pb)
  }

  cli::cli_progress_done(id = pb)
  log_msg("{length(downloaded)}/{nrow(scenes)} scènes S2 téléchargées dans {s2_dir}",
          level = "success")

  extract_s2_bands(s2_dir, aoi)

  s2_dir
}

#' Téléchargement de tous les produits S1 pour une AOI et une année
#' @param aoi sf object
#' @param year Année
#' @param output_dir Répertoire de sortie
#' @param token Token CDSE
#' @param max_scenes Nombre max de scènes
#' @return Chemin vers le répertoire des données
download_s1_for_aoi <- function(aoi, year = 2021,
                                 output_dir = file.path(RAW_DIR, "sentinel1"),
                                 token = NULL, max_scenes = NULL) {
  cli::cli_h2("Téléchargement Sentinel-1 GRD")

  if (is.null(token)) token <- cdse_get_token()
  if (is.null(token)) return(NULL)

  start_date <- paste0(year, "-01-01")
  end_date   <- paste0(year, "-12-31")
  scenes <- search_sentinel1(aoi, start_date, end_date, orbit_direction = "DESCENDING")
  if (is.null(scenes)) return(NULL)

  # Sélection : 1 scène tous les ~12 jours
  scenes <- select_best_scenes_s1(scenes, interval_days = 12)
  log_msg("  {nrow(scenes)} scènes S1 sélectionnées")

  if (!is.null(max_scenes)) {
    scenes <- scenes[1:min(max_scenes, nrow(scenes)), ]
  }

  if (interactive()) {
    confirm <- readline(glue::glue(
      "Télécharger {nrow(scenes)} scènes S1 ? (o/n) : "
    ))
    if (!tolower(confirm) %in% c("o", "oui", "y", "yes")) {
      log_msg("Téléchargement S1 annulé", level = "warning")
      return(NULL)
    }
  }

  s1_dir <- file.path(output_dir, paste0("S1_GRD_", year))
  dir.create(s1_dir, showWarnings = FALSE, recursive = TRUE)

  pb <- cli::cli_progress_bar("Téléchargement S1", total = nrow(scenes))

  downloaded <- c()
  for (i in seq_len(nrow(scenes))) {
    scene <- scenes[i, ]
    log_msg("  [{i}/{nrow(scenes)}] {scene$name} ({scene$date})")

    zip_path <- download_product(scene$id, scene$name, token, s1_dir)
    if (!is.null(zip_path)) {
      if (!dir.exists(file.path(s1_dir, tools::file_path_sans_ext(basename(zip_path))))) {
        unzip(zip_path, exdir = s1_dir)
      }
      downloaded <- c(downloaded, zip_path)
    }

    cli::cli_progress_update(id = pb)
  }

  cli::cli_progress_done(id = pb)
  log_msg("{length(downloaded)}/{nrow(scenes)} scènes S1 téléchargées dans {s1_dir}",
          level = "success")

  s1_dir
}

# ==============================================================================
# 4. PIPELINE COMPLET DE TÉLÉCHARGEMENT S2 + S1
# ==============================================================================

#' Téléchargement automatique de toutes les données satellite pour une AOI
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
  cli::cli_h1("Téléchargement des données satellite")
  t_start <- Sys.time()

  aoi <- sf::st_read(aoi_path, quiet = TRUE)
  if (nrow(aoi) > 1) aoi <- sf::st_union(aoi) |> sf::st_as_sf()

  aoi_area_ha <- as.numeric(sf::st_area(sf::st_transform(aoi, 2154))) / 10000
  log_msg("AOI : {round(aoi_area_ha, 1)} ha")

  token <- cdse_get_token()
  if (is.null(token)) {
    cli::cli_alert_danger("Impossible de s'authentifier — abandon")
    return(NULL)
  }

  results <- list(s2_dir = NULL, s1_dir = NULL)

  if (download_s2) {
    results$s2_dir <- download_s2_for_aoi(
      aoi, year = year, max_cloud = max_cloud,
      output_dir = output_dir, token = token
    )
  }

  if (download_s1) {
    results$s1_dir <- download_s1_for_aoi(
      aoi, year = year,
      output_dir = output_dir, token = token
    )
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  cli::cli_h2("Téléchargement terminé ({round(t_elapsed, 1)} min)")

  if (!is.null(results$s2_dir)) {
    n_s2 <- length(list.files(results$s2_dir, pattern = "\\.zip$", recursive = TRUE))
    cli::cli_alert_success("Sentinel-2 : {n_s2} produits → {results$s2_dir}")
  }
  if (!is.null(results$s1_dir)) {
    n_s1 <- length(list.files(results$s1_dir, pattern = "\\.zip$", recursive = TRUE))
    cli::cli_alert_success("Sentinel-1 : {n_s1} produits → {results$s1_dir}")
  }

  invisible(results)
}

# ==============================================================================
# 5. FONCTIONS UTILITAIRES
# ==============================================================================

#' Sélection des meilleures scènes S2 (1 par intervalle, moins de nuages)
#' @param scenes data.frame de scènes
#' @param interval_days Intervalle cible entre scènes (jours)
#' @return data.frame des scènes sélectionnées
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

#' Extraction et réorganisation des bandes S2 depuis les archives SAFE
#' @param s2_dir Répertoire contenant les produits S2
#' @param aoi sf object pour le crop
extract_s2_bands <- function(s2_dir, aoi) {
  log_msg("Extraction des bandes S2...")

  bands_dir <- file.path(s2_dir, "bands")
  dir.create(bands_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_proj <- sf::st_transform(aoi, 2154)
  aoi_vect <- terra::vect(aoi_proj)
  aoi_ext  <- terra::ext(terra::buffer(aoi_vect, 500))

  all_bands <- list.files(s2_dir, pattern = "(B02|B03|B04|B05|B06|B07|B08|B8A|B11|B12|SCL).*\\.(jp2|tif)$",
                          recursive = TRUE, full.names = TRUE)

  if (length(all_bands) == 0) {
    log_msg("  Aucune bande trouvée à extraire", level = "warning")
    return(invisible(NULL))
  }

  log_msg("  {length(all_bands)} fichiers de bandes trouvés")

  pb <- cli::cli_progress_bar("Extraction bandes", total = length(all_bands))

  for (f in all_bands) {
    band_name <- stringr::str_extract(basename(f), "B\\d{2}|B8A|SCL")
    date_str  <- stringr::str_extract(f, "\\d{8}T\\d{6}")
    if (is.na(date_str)) date_str <- stringr::str_extract(basename(f), "\\d{8}")

    if (is.na(band_name) || is.na(date_str)) {
      cli::cli_progress_update(id = pb)
      next
    }

    date_short <- substr(date_str, 1, 8)
    out_name <- paste0(date_short, "_", band_name, ".tif")
    out_path <- file.path(bands_dir, out_name)

    if (!file.exists(out_path)) {
      tryCatch({
        r <- terra::rast(f)
        if (!terra::same.crs(r, terra::crs("EPSG:2154"))) {
          r <- terra::project(r, "EPSG:2154", method = "bilinear")
        }
        r <- terra::crop(r, aoi_ext)
        terra::writeRaster(r, out_path, overwrite = TRUE)
      }, error = function(e) {
      })
    }

    cli::cli_progress_update(id = pb)
  }

  cli::cli_progress_done(id = pb)

  n_extracted <- length(list.files(bands_dir, pattern = "\\.tif$"))
  log_msg("  {n_extracted} bandes extraites et croppées dans {bands_dir}", level = "success")

  invisible(bands_dir)
}

#' Prétraitement des données Sentinel-1 GRD
#' @param s1_dir Répertoire contenant les produits S1
#' @param aoi sf object pour le crop
#' @param resolution Résolution cible en mètres
#' @return Chemin vers les bandes S1 traitées
preprocess_s1 <- function(s1_dir, aoi, resolution = 10) {
  log_msg("Prétraitement Sentinel-1 GRD...")

  bands_dir <- file.path(s1_dir, "bands")
  dir.create(bands_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_proj <- sf::st_transform(aoi, 2154)
  aoi_vect <- terra::vect(aoi_proj)
  aoi_ext  <- terra::ext(terra::buffer(aoi_vect, 500))

  tiff_files <- list.files(s1_dir, pattern = "(vv|vh).*\\.tiff?$",
                           recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  if (length(tiff_files) == 0) {
    log_msg("  Aucun fichier S1 trouvé", level = "warning")
    return(invisible(NULL))
  }

  log_msg("  {length(tiff_files)} fichiers S1 trouvés")

  pb <- cli::cli_progress_bar("Traitement S1", total = length(tiff_files))

  for (f in tiff_files) {
    pol <- toupper(stringr::str_extract(basename(f), "(?i)(vv|vh)"))
    date_str <- stringr::str_extract(f, "\\d{8}T\\d{6}")
    if (is.na(date_str)) date_str <- stringr::str_extract(basename(f), "\\d{8}")

    if (is.na(pol) || is.na(date_str)) {
      cli::cli_progress_update(id = pb)
      next
    }

    date_short <- substr(date_str, 1, 8)
    out_name <- paste0(date_short, "_", pol, "_sigma0_db.tif")
    out_path <- file.path(bands_dir, out_name)

    if (!file.exists(out_path)) {
      tryCatch({
        r <- terra::rast(f)

        if (!terra::same.crs(r, terra::crs("EPSG:2154"))) {
          r <- terra::project(r, "EPSG:2154", method = "bilinear",
                              res = resolution)
        }

        r <- terra::crop(r, aoi_ext)

        # Calibration : DN → sigma0 en dB
        vals <- terra::values(r)
        vals[vals <= 0] <- NA
        vals <- 10 * log10(vals)
        terra::values(r) <- vals

        # Filtre de Lee (speckle)
        r <- terra::focal(r, w = matrix(1, 3, 3), fun = "mean", na.rm = TRUE)

        r <- terra::mask(r, aoi_vect)

        terra::writeRaster(r, out_path, overwrite = TRUE)
      }, error = function(e) {
      })
    }

    cli::cli_progress_update(id = pb)
  }

  cli::cli_progress_done(id = pb)

  n_processed <- length(list.files(bands_dir, pattern = "\\.tif$"))
  log_msg("  {n_processed} bandes S1 traitées dans {bands_dir}", level = "success")

  invisible(bands_dir)
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
    log_msg("  Aucun fichier S1 traité trouvé", level = "warning")
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
#' @param vv Valeur VV en dB
#' @param vh Valeur VH en dB
#' @return Liste d'indices radar
calc_radar_indices <- function(vv, vh) {
  list(
    VV_VH_ratio = vv - vh,
    RVI = 4 * 10^(vh/10) / (10^(vv/10) + 10^(vh/10)),
    RFDI = (10^(vv/10) - 10^(vh/10)) /
           (10^(vv/10) + 10^(vh/10))
  )
}

#' Statistiques temporelles Sentinel-1 pour un pixel
#' @param vv_ts Série temporelle VV (dB)
#' @param vh_ts Série temporelle VH (dB)
#' @return Vecteur nommé de features S1
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

cli::cli_alert_success("Module de téléchargement satellite chargé (rstac)")
cli::cli_text("Utilisation :")
cli::cli_text('  {.code download_satellite_data("aoi.gpkg", year = 2023)}')
cli::cli_text('  {.code download_s2_for_aoi(aoi, year = 2023)}')
cli::cli_text('  {.code download_s1_for_aoi(aoi, year = 2023)}')
