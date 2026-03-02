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
# 0. MASQUE FORESTIER (OSO + NDVI)
# ==============================================================================

#' Télécharger les données d'occupation du sol pour une AOI
#'
#' Stratégie multi-source (du plus léger au plus lourd) :
#'
#'   1. **WFS Géoplateforme IGN** (prioritaire, quelques Ko) :
#'      BD Forêt V2 : couche `LANDCOVER.FORESTINVENTORY.V2:formation_vegetale`
#'      Accès libre, données vectorielles découpées à la bbox, rasterisées sur place.
#'      Ref : https://geoservices.ign.fr/services-web-experts-ocsge
#'
#'   2. **OCS GE par département** (GPKG, quelques dizaines de Mo) :
#'      Téléchargement via l'API Géoplateforme, archives .7z par département.
#'      Nomenclature couverture : CS2.1.1.1 = Feuillus, CS2.1.1.2 = Conifères,
#'      CS2.1.1.3 = Mixte. Ref : https://geoservices.ign.fr/ocsge#telechargement
#'
#'   3. **OSO raster CESBIO** (~6 Go, cache global Recherche Data Gouv) :
#'      Raster 10 m France entière, classes 16 = Feuillus / 17 = Conifères.
#'      Ref : https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/UZ2NJ7
#'
#'   4. **Fichier local** : recherche dans data/raw/, data/, cache global.
#'
#' @param aoi sf object — zone d'intérêt
#' @param year Année (pour la correspondance millésime OCS GE / OSO)
#' @param output_dir Répertoire de sortie (cache projet)
#' @param resolution Résolution cible pour la rastérisation (mètres, défaut 10)
#' @return Chemin vers le raster d'occupation du sol découpé à l'AOI (GeoTIFF)
#' @export
download_oso <- function(aoi, year = FOREST_MASK_PARAMS$oso_year,
                          output_dir = file.path(RAW_DIR, "oso"),
                          resolution = 10) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  oso_file <- file.path(output_dir, paste0("oso_", year, ".tif"))

  # Retourner si déjà découpé pour cette AOI
  if (file.exists(oso_file)) {
    log_msg("  OSO {year} déjà présent : {oso_file}", level = "info")
    return(oso_file)
  }

  log_msg("Acquisition du masque d'occupation du sol...")

  # L'AOI en Lambert-93 et en WGS84
  aoi_proj <- sf::st_transform(aoi, 2154)
  aoi_wgs84 <- sf::st_transform(aoi, 4326)
  bbox_wgs84 <- sf::st_bbox(aoi_wgs84)

  # =====================================================================
  # Stratégie 1 : WFS Géoplateforme IGN (BD Forêt V2 — vecteur, léger)
  # =====================================================================
  # La BD Forêt V2 contient les formations végétales avec le type de forêt.
  # Accès libre via WFS, données vectorielles, quelques Ko pour une AOI.
  # Même approche que le package nemeton (download_ign_bdforet).

  log_msg("  Tentative WFS Géoplateforme (BD Forêt V2)...")
  wfs_ok <- tryCatch({
    wfs_url <- "https://data.geopf.fr/wfs/ows"
    typename <- "LANDCOVER.FORESTINVENTORY.V2:formation_vegetale"

    # Construire la requête WFS GetFeature avec bbox
    bbox_str <- paste(
      round(bbox_wgs84["ymin"], 6), round(bbox_wgs84["xmin"], 6),
      round(bbox_wgs84["ymax"], 6), round(bbox_wgs84["xmax"], 6),
      sep = ","
    )

    # --- Pagination WFS : récupérer TOUS les polygones ---
    # Le serveur Géoplateforme limite souvent à 1000-5000 features par
    # requête. On pagine avec STARTINDEX + COUNT jusqu'à épuisement.
    page_size <- 5000L
    start_idx <- 0L
    all_features <- list()

    repeat {
      wfs_request <- paste0(
        wfs_url,
        "?SERVICE=WFS",
        "&VERSION=2.0.0",
        "&REQUEST=GetFeature",
        "&TYPENAMES=", typename,
        "&BBOX=", bbox_str, ",EPSG:4326",
        "&OUTPUTFORMAT=application/json",
        "&COUNT=", page_size,
        "&STARTINDEX=", start_idx
      )

      tmp_json <- tempfile(fileext = ".json")
      resp <- httr2::request(wfs_request) |>
        httr2::req_timeout(120) |>
        httr2::req_perform()

      httr2::resp_body_raw(resp) |> writeBin(tmp_json)
      page_sf <- sf::st_read(tmp_json, quiet = TRUE)
      unlink(tmp_json)

      if (nrow(page_sf) == 0) break

      all_features <- c(all_features, list(page_sf))
      log_msg("  WFS page : {start_idx + nrow(page_sf)} polygones récupérés...")

      # Si on a reçu moins que page_size, c'est la dernière page
      if (nrow(page_sf) < page_size) break
      start_idx <- start_idx + page_size
    }

    if (length(all_features) == 0) {
      log_msg("  Aucune formation végétale trouvée via WFS", level = "warning")
      FALSE
    } else {
      bdforet <- do.call(rbind, all_features)
      log_msg("  BD Forêt V2 : {nrow(bdforet)} polygones récupérés (total)", level = "success")

      # Reprojeter en Lambert-93
      bdforet <- sf::st_transform(bdforet, 2154)

      # Créer un raster template aligné sur l'AOI
      aoi_ext <- terra::ext(terra::vect(aoi_proj))
      template <- terra::rast(aoi_ext, resolution = resolution, crs = "EPSG:2154")

      # Rastériser : toutes les formations végétales = 1 (forêt)
      # La BD Forêt V2 ne contient QUE des zones boisées
      bdforet_vect <- terra::vect(bdforet)
      r_forest <- terra::rasterize(bdforet_vect, template, field = 1, background = 0)
      names(r_forest) <- "forest"

      # Masquer à l'AOI
      aoi_vect <- terra::vect(aoi_proj)
      r_forest <- terra::mask(r_forest, aoi_vect)

      terra::writeRaster(r_forest, oso_file, datatype = "INT1U", overwrite = TRUE)
      log_msg("  Masque forêt (BD Forêt V2) sauvegardé : {oso_file}", level = "success")
      TRUE
    }
  }, error = function(e) {
    log_msg("  WFS Géoplateforme indisponible : {e$message}", level = "warning")
    FALSE
  })

  if (wfs_ok) return(oso_file)

  # =====================================================================
  # Stratégie 2 : OCS GE par département (GPKG via API Géoplateforme)
  # =====================================================================
  # L'OCS GE (Occupation du Sol à Grande Échelle, IGN) est disponible en
  # téléchargement par département au format GeoPackage (.7z).
  # Nomenclature couverture du sol (CS) :
  #   CS2.1.1.1 = Peuplements de feuillus
  #   CS2.1.1.2 = Peuplements de conifères
  #   CS2.1.1.3 = Peuplements mixtes
  # Ref : https://geoservices.ign.fr/ocsge#telechargement

  log_msg("  Tentative OCS GE par département (Géoplateforme)...")
  ocsge_ok <- tryCatch({
    # --- Déterminer TOUS les départements couverts par l'AOI ---
    # Échantillonner le centroïde + les coins de la bbox pour détecter
    # les AOIs chevauchant plusieurs départements.
    bbox_w <- bbox_wgs84["xmin"]; bbox_e <- bbox_wgs84["xmax"]
    bbox_s <- bbox_wgs84["ymin"]; bbox_n <- bbox_wgs84["ymax"]
    sample_pts <- data.frame(
      lon = c(
        (bbox_w + bbox_e) / 2,  # centroïde
        bbox_w, bbox_e, bbox_w, bbox_e  # 4 coins
      ),
      lat = c(
        (bbox_s + bbox_n) / 2,
        bbox_s, bbox_s, bbox_n, bbox_n
      )
    )

    # Fonction interne : code département depuis coordonnées WGS84
    .get_dept_code <- function(lon, lat) {
      tryCatch({
        rev_url <- paste0(
          "https://data.geopf.fr/geocodage/reverse?lon=", lon,
          "&lat=", lat, "&type=municipality&limit=1"
        )
        rev_resp <- httr2::request(rev_url) |>
          httr2::req_timeout(10) |>
          httr2::req_perform()
        rev_json <- httr2::resp_body_json(rev_resp)
        if (length(rev_json$features) == 0) return(NULL)
        props <- rev_json$features[[1]]$properties
        code_insee <- props$citycode %||% props$postcode %||% ""
        if (nchar(code_insee) < 5) return(NULL)
        if (startsWith(code_insee, "97")) {
          substr(code_insee, 1, 3)
        } else if (startsWith(code_insee, "20")) {
          cc <- as.integer(substr(code_insee, 3, 5))
          if (!is.na(cc) && cc >= 1 && cc <= 360) "2A" else "2B"
        } else {
          substr(code_insee, 1, 2)
        }
      }, error = function(e) NULL)
    }

    # Récupérer les départements uniques
    dept_codes <- unique(Filter(Negate(is.null), lapply(
      seq_len(nrow(sample_pts)),
      function(i) .get_dept_code(sample_pts$lon[i], sample_pts$lat[i])
    )))
    dept_codes <- unique(unlist(dept_codes))

    if (length(dept_codes) == 0) {
      log_msg("  Impossible de déterminer le(s) département(s)", level = "warning")
      stop("département inconnu")
    }

    log_msg("  Département(s) détecté(s) : {paste(dept_codes, collapse = ', ')}")

    # Fonction interne : formater le code département pour l'API
    .format_dept_zone <- function(dc) {
      if (grepl("^[0-9]+$", dc)) {
        paste0("D", formatC(as.integer(dc), width = 3, flag = "0"))
      } else {
        paste0("D0", dc)  # D02A, D02B
      }
    }

    # --- Télécharger et fusionner les polygones forestiers de chaque département ---
    ocsge_cache <- file.path(
      rappdirs::user_cache_dir("treesatnemeton"), "ocsge"
    )
    dir.create(ocsge_cache, showWarnings = FALSE, recursive = TRUE)

    all_forest_polys <- list()

    for (dept_code in dept_codes) {
      dept_zone <- .format_dept_zone(dept_code)
      log_msg("  Traitement OCS GE département {dept_code} (zone {dept_zone})...")

      # Interroger l'API Géoplateforme pour trouver les fichiers OCS GE
      dataset_name <- NULL
      for (fmt in c("GPKG", "SHP")) {
        api_url <- paste0(
          "https://data.geopf.fr/telechargement/resource/OCSGE",
          "?format=", fmt, "&zone=", dept_zone, "&page=1&limit=50"
        )
        api_resp <- tryCatch(
          httr2::request(api_url) |>
            httr2::req_timeout(30) |>
            httr2::req_perform(),
          error = function(e) NULL
        )
        if (is.null(api_resp)) next
        api_json <- httr2::resp_body_json(api_resp)
        entries <- api_json$datasets %||% api_json$entries %||% api_json
        for (entry in entries) {
          nm <- entry$name %||% entry$title %||% ""
          if (grepl(fmt, nm) && !grepl("DIFF", nm) && !grepl("ARTIF", nm)) {
            dataset_name <- nm
            break
          }
        }
        if (!is.null(dataset_name)) break
      }

      if (is.null(dataset_name)) {
        log_msg("  OCS GE non disponible pour {dept_code}", level = "warning")
        next
      }

      log_msg("  OCS GE trouvé : {dataset_name}")

      # Télécharger l'archive .7z
      archive_file <- file.path(ocsge_cache, paste0(dataset_name, ".7z"))

      if (!file.exists(archive_file)) {
        dl_url <- paste0(
          "https://data.geopf.fr/telechargement/download/OCSGE/",
          dataset_name, "/", dataset_name, ".7z"
        )
        log_msg("  Téléchargement OCS GE {dept_code}...")
        utils::download.file(dl_url, archive_file, mode = "wb", quiet = TRUE)
      }

      # Extraire l'archive .7z
      extract_dir <- file.path(ocsge_cache, dataset_name)
      if (!dir.exists(extract_dir)) {
        log_msg("  Extraction de l'archive .7z...")
        if (requireNamespace("archive", quietly = TRUE)) {
          archive::archive_extract(archive_file, dir = extract_dir)
        } else {
          sys_7z <- Sys.which("7z")
          if (nchar(sys_7z) == 0) sys_7z <- Sys.which("7za")
          if (nchar(sys_7z) == 0) {
            log_msg("  Package 'archive' ou commande '7z' requis pour OCS GE .7z",
                    level = "warning")
            next
          }
          dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
          system2(sys_7z, args = c("x", "-y", paste0("-o", extract_dir),
                                    archive_file), stdout = FALSE, stderr = FALSE)
        }
      }

      # Trouver le fichier GPKG ou SHP
      gpkg_files <- list.files(extract_dir, pattern = "\\.gpkg$",
                                recursive = TRUE, full.names = TRUE)
      shp_files <- list.files(extract_dir, pattern = "\\.shp$",
                               recursive = TRUE, full.names = TRUE)

      ocsge_sf <- NULL
      if (length(gpkg_files) > 0) {
        lyrs <- sf::st_layers(gpkg_files[1])$name
        couv_lyr <- lyrs[grepl("couverture", lyrs, ignore.case = TRUE)]
        if (length(couv_lyr) == 0) couv_lyr <- lyrs[1]
        ocsge_sf <- sf::st_read(gpkg_files[1], layer = couv_lyr[1], quiet = TRUE)
      } else if (length(shp_files) > 0) {
        couv_shp <- shp_files[grepl("couverture", shp_files, ignore.case = TRUE)]
        if (length(couv_shp) == 0) couv_shp <- shp_files[1]
        ocsge_sf <- sf::st_read(couv_shp[1], quiet = TRUE)
      }

      if (is.null(ocsge_sf) || nrow(ocsge_sf) == 0) {
        log_msg("  Aucune donnée couverture pour {dept_code}", level = "warning")
        next
      }

      log_msg("  OCS GE {dept_code} : {nrow(ocsge_sf)} polygones chargés")

      # Filtrer les formations arborées (CS2.1.1.*)
      cs_col <- intersect(names(ocsge_sf), c("code_cs", "CODE_CS", "couverture",
                                              "COUVERTURE", "code_couv"))
      if (length(cs_col) == 0) {
        for (cn in names(ocsge_sf)) {
          if (is.character(ocsge_sf[[cn]])) {
            sample_vals <- head(stats::na.omit(ocsge_sf[[cn]]), 20)
            if (any(grepl("^CS", sample_vals))) {
              cs_col <- cn
              break
            }
          }
        }
      }
      if (length(cs_col) == 0) {
        log_msg("  Colonne code_cs non trouvée pour {dept_code}", level = "warning")
        next
      }
      cs_col <- cs_col[1]

      forest_pattern <- "^CS2\\.1\\.1"
      forest_polys <- ocsge_sf[grepl(forest_pattern, ocsge_sf[[cs_col]]), ]

      if (nrow(forest_polys) > 0) {
        forest_polys <- sf::st_transform(forest_polys, 2154)
        all_forest_polys <- c(all_forest_polys, list(forest_polys))
        log_msg("  {nrow(forest_polys)} polygones forestiers pour {dept_code}")
      }
    }  # fin boucle départements

    if (length(all_forest_polys) == 0) {
      log_msg("  Aucune formation arborée (CS2.1.1.*) dans la zone", level = "warning")
      stop("pas de forêt OCS GE")
    }

    # Fusionner les polygones de tous les départements
    if (length(all_forest_polys) > 1) {
      # Uniformiser les colonnes avant rbind
      common_cols <- Reduce(intersect, lapply(all_forest_polys, names))
      all_forest_polys <- lapply(all_forest_polys, function(x) x[, common_cols])
      forest_merged <- do.call(rbind, all_forest_polys)
    } else {
      forest_merged <- all_forest_polys[[1]]
    }

    log_msg("  Total : {nrow(forest_merged)} polygones forestiers ({length(dept_codes)} département(s))")

    # Découper à l'AOI
    forest_merged <- sf::st_intersection(forest_merged, sf::st_geometry(aoi_proj))

    # Rastériser
    aoi_ext <- terra::ext(terra::vect(aoi_proj))
    template <- terra::rast(aoi_ext, resolution = resolution, crs = "EPSG:2154")
    forest_vect <- terra::vect(forest_merged)
    r_forest <- terra::rasterize(forest_vect, template, field = 1, background = 0)
    names(r_forest) <- "forest"

    aoi_vect <- terra::vect(aoi_proj)
    r_forest <- terra::mask(r_forest, aoi_vect)

    terra::writeRaster(r_forest, oso_file, datatype = "INT1U", overwrite = TRUE)
    log_msg("  Masque forêt (OCS GE {paste(dept_codes, collapse='/')}) sauvegardé : {oso_file}",
            level = "success")
    TRUE
  }, error = function(e) {
    log_msg("  OCS GE par département échoué : {e$message}", level = "warning")
    FALSE
  })

  if (ocsge_ok) return(oso_file)

  # =====================================================================
  # Stratégie 3 : OSO raster CESBIO (Recherche Data Gouv, ~6 Go)
  # =====================================================================
  log_msg("  Tentative OSO raster (Recherche Data Gouv)...")

  # Cache global : raster France entière (partagé entre projets)
  global_cache <- file.path(
    rappdirs::user_cache_dir("treesatnemeton"), "oso"
  )
  dir.create(global_cache, showWarnings = FALSE, recursive = TRUE)
  oso_global <- file.path(global_cache, "oso.tif")

  if (!file.exists(oso_global)) {
    # Téléchargement depuis Recherche Data Gouv
    # Source : https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/UZ2NJ7
    oso_url <- "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/8M1AN1"

    oso_tar <- file.path(global_cache, "OSO_RASTER.tar.gz")
    oso_tar_files <- list.files(global_cache, pattern = "^OSO_.*\\.tar\\.gz$",
                                 full.names = TRUE)
    if (length(oso_tar_files) > 0) oso_tar <- oso_tar_files[1]

    need_download <- TRUE
    oso_expected_size <- 5e9

    if (file.exists(oso_tar) && file.info(oso_tar)$size >= oso_expected_size) {
      log_msg("  Archive OSO existante ({round(file.info(oso_tar)$size/1e9, 1)} Go)",
              level = "info")
      need_download <- FALSE
    }

    if (need_download) {
      log_msg("  Téléchargement OSO (~6 Go, patience...)...")
      log_msg("  Cache global : {global_cache}", level = "info")

      download_ok <- tryCatch({
        old_timeout <- getOption("timeout")
        options(timeout = 3600)
        on.exit(options(timeout = old_timeout), add = TRUE)

        if (requireNamespace("curl", quietly = TRUE)) {
          h <- curl::new_handle()
          curl::handle_setopt(h,
            followlocation = TRUE, timeout = 3600,
            low_speed_limit = 1000, low_speed_time = 300
          )
          con <- curl::curl(oso_url, handle = h, open = "rb")
          on.exit(try(close(con), silent = TRUE), add = TRUE)
          out <- file(oso_tar, open = "wb")
          on.exit(try(close(out), silent = TRUE), add = TRUE)

          downloaded <- 0
          chunk_size <- 1024 * 1024
          last_pct <- -1
          while (TRUE) {
            buf <- readBin(con, raw(), n = chunk_size)
            if (length(buf) == 0) break
            writeBin(buf, out)
            downloaded <- downloaded + length(buf)
            pct <- min(floor(downloaded / oso_expected_size * 100), 99)
            if (pct %% 10 == 0 && pct != last_pct) {
              log_msg("  Téléchargement : {pct}% ({round(downloaded/1e9, 1)} Go)")
              last_pct <- pct
            }
          }
          close(out); close(con)
          log_msg("  OSO téléchargé : {round(downloaded/1e9, 1)} Go", level = "success")
          TRUE
        } else {
          utils::download.file(oso_url, oso_tar, mode = "wb", quiet = FALSE)
          TRUE
        }
      }, error = function(e) {
        log_msg("  Échec téléchargement OSO : {e$message}", level = "warning")
        FALSE
      })
    }

    # Extraire l'archive
    if (file.exists(oso_tar) && !file.exists(oso_global)) {
      log_msg("  Extraction de l'archive OSO...")
      utils::untar(oso_tar, exdir = global_cache)
      tif_files <- list.files(global_cache, pattern = "^OCS_.*\\.tif$",
                               recursive = TRUE, full.names = TRUE)
      if (length(tif_files) > 0) {
        file.copy(tif_files[1], oso_global, overwrite = TRUE)
        log_msg("  OSO extrait : {oso_global}", level = "success")
        oso_dirs <- list.dirs(global_cache, full.names = TRUE, recursive = FALSE)
        oso_dirs <- oso_dirs[grepl("^OSO_", basename(oso_dirs))]
        if (length(oso_dirs) > 0) unlink(oso_dirs, recursive = TRUE)
      }
    }
  }

  # =====================================================================
  # Stratégie 4 : Fichier local existant
  # =====================================================================
  if (!file.exists(oso_global)) {
    search_dirs <- unique(c(RAW_DIR, DATA_DIR, output_dir, global_cache))
    search_dirs <- search_dirs[dir.exists(search_dirs)]
    existing <- list.files(
      search_dirs,
      pattern = "(OSO|oso|OCS).*\\.(tif|TIF)$",
      recursive = TRUE, full.names = TRUE
    )
    if (length(existing) > 0) {
      log_msg("  OSO trouvé localement : {existing[1]}", level = "info")
      oso_global <- existing[1]
    }
  }

  if (!file.exists(oso_global)) {
    log_msg("  Aucune donnée d'occupation du sol disponible.", level = "warning")
    log_msg("  Masque NDVI seul sera utilisé.", level = "warning")
    return(NULL)
  }

  # Découper le raster OSO à l'AOI
  # IMPORTANT : on découpe à la bbox élargie (+ 500 m de marge) et NON au

  # polygone exact, pour que le rééchantillonnage sur la grille S2 dans
  # build_oso_forest_mask() ne produise pas de NA en bordure.
  log_msg("  Découpe OSO à l'AOI (bbox + marge)...")
  tryCatch({
    r <- terra::rast(oso_global)
    aoi_vect <- terra::vect(aoi_proj)
    if (!terra::same.crs(r, aoi_vect)) {
      aoi_vect <- terra::project(aoi_vect, terra::crs(r))
    }
    # Étendre la bbox de 500 m pour couvrir les pixels de bordure du cube S2
    crop_ext <- terra::ext(aoi_vect)
    crop_ext[1] <- crop_ext[1] - 500  # xmin
    crop_ext[2] <- crop_ext[2] + 500  # xmax
    crop_ext[3] <- crop_ext[3] - 500  # ymin
    crop_ext[4] <- crop_ext[4] + 500  # ymax
    r_crop <- terra::crop(r, crop_ext)
    terra::writeRaster(r_crop, oso_file, datatype = "INT1U", overwrite = TRUE)
    log_msg("  OSO découpé sauvegardé : {oso_file}", level = "success")
    return(oso_file)
  }, error = function(e) {
    log_msg("  Erreur découpe OSO : {e$message}", level = "warning")
    return(NULL)
  })
}

#' Construire un masque forêt à partir du raster d'occupation du sol (binaire)
#'
#' Convertit le raster d'occupation du sol en masque binaire forêt (1) / non-forêt (0),
#' rééchantillonné sur la grille du cube Sentinel-2.
#'
#' Gère trois cas :
#'   - BD Forêt V2 / OCS GE (via WFS ou téléchargement) : déjà binaire (1 = forêt, 0 = non-forêt)
#'   - OSO CESBIO : classes 16 = Feuillus, 17 = Conifères → binaire
#'
#' @param oso_path Chemin vers le raster d'occupation du sol
#' @param template SpatRaster — grille cible (même résolution/emprise que le cube S2)
#' @param forest_classes Codes OSO considérés comme forestiers (ignoré si BD Forêt)
#' @return SpatRaster binaire aligné sur le template (1 = forêt, 0 = non-forêt)
#' @export
build_oso_forest_mask <- function(oso_path, template,
                                   forest_classes = FOREST_MASK_PARAMS$oso_forest_classes) {
  r_oso <- terra::rast(oso_path)

  # Vérifier la couverture spatiale OSO vs template S2
  oso_ext <- terra::ext(r_oso)
  tpl_ext <- terra::ext(template)
  if (!terra::same.crs(r_oso, template)) {
    # Reprojeter si CRS différents
    r_oso <- terra::project(r_oso, template, method = "near")
    oso_ext <- terra::ext(r_oso)
  }

  # Diagnostic de couverture
  n_total <- terra::ncell(template)
  covers_fully <- (oso_ext[1] <= tpl_ext[1]) && (oso_ext[2] >= tpl_ext[2]) &&
                  (oso_ext[3] <= tpl_ext[3]) && (oso_ext[4] >= tpl_ext[4])

  if (!covers_fully) {
    log_msg("  Attention : le raster OSO ne couvre pas entièrement le template S2", level = "warning")
    log_msg("  OSO ext : {round(oso_ext[1])}, {round(oso_ext[2])}, {round(oso_ext[3])}, {round(oso_ext[4])}")
    log_msg("  S2  ext : {round(tpl_ext[1])}, {round(tpl_ext[2])}, {round(tpl_ext[3])}, {round(tpl_ext[4])}")
    # Étendre le raster OSO à l'emprise du template (NA pour les pixels manquants)
    r_oso <- terra::extend(r_oso, template)
  }

  # Rééchantillonner sur la grille S2 (nearest neighbor pour catégoriel)
  r_oso <- terra::resample(r_oso, template, method = "near")

  oso_vals <- terra::values(r_oso)[, 1]
  unique_vals <- unique(stats::na.omit(oso_vals))

  # Détecter si le raster est déjà binaire (BD Forêt WFS : valeurs 0/1 uniquement)
  is_binary <- length(unique_vals) <= 2 && all(unique_vals %in% c(0, 1))

  if (is_binary) {
    # BD Forêt V2 : déjà 1 = forêt, 0 = non-forêt
    mask_vals <- as.integer(oso_vals == 1L)
  } else {
    # OSO CESBIO : filtrer par classes forestières (16 = Feuillus, 17 = Conifères)
    mask_vals <- as.integer(oso_vals %in% forest_classes)
  }
  mask_vals[is.na(oso_vals)] <- 0L

  # Diagnostic couverture
  n_na <- sum(is.na(oso_vals))
  if (n_na > 0) {
    pct_na <- round(n_na / n_total * 100, 1)
    log_msg("  OSO : {n_na}/{n_total} pixels sans donnée ({pct_na}%) → traités comme non-forêt",
            level = if (pct_na > 5) "warning" else "info")
  }

  r_mask <- terra::rast(template)
  terra::values(r_mask) <- mask_vals
  names(r_mask) <- "forest_oso"
  r_mask
}

#' Construire un masque NDVI (végétation active)
#'
#' Identifie les pixels avec un NDVI max annuel supérieur au seuil.
#' Détecte la végétation active même si la carte OSO est obsolète
#' (jeune régénération, forêt récente hors carte OSO).
#'
#' @param cube_arrays Liste des matrices par bande (n_pixels × n_dates)
#' @param template SpatRaster — grille cible
#' @param threshold Seuil NDVI max annuel (défaut 0.4)
#' @return SpatRaster binaire (1 = NDVI max >= seuil, 0 sinon)
#' @export
build_ndvi_mask <- function(cube_arrays, template,
                             threshold = FOREST_MASK_PARAMS$ndvi_min_threshold) {
  # Calcul du NDVI pour toutes les dates
  b08 <- cube_arrays[["B08"]]
  b04 <- cube_arrays[["B04"]]

  # NDVI = (NIR - Red) / (NIR + Red)
  ndvi <- (b08 - b04) / (b08 + b04 + 1e-10)

  # NDVI max annuel par pixel (max sur toutes les dates, en ignorant les NA)
  ndvi_max <- apply(ndvi, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })

  # Masque binaire : 1 si NDVI max >= seuil
  mask_vals <- as.integer(!is.na(ndvi_max) & ndvi_max >= threshold)

  r_mask <- terra::rast(template)
  terra::values(r_mask) <- mask_vals
  names(r_mask) <- "forest_ndvi"
  r_mask
}

#' Construire le masque forestier combiné (OSO + NDVI)
#'
#' Combine les deux sources pour créer un masque robuste :
#' - "union" : pixel forestier si OSO forêt OU NDVI > seuil
#'   (conserve les coupes rases avec NDVI résiduel et les jeunes plantations)
#' - "intersection" : pixel forestier si OSO forêt ET NDVI > seuil
#'   (plus strict, exclut les coupes rases dans les zones OSO « forêt »)
#'
#' @param oso_mask SpatRaster binaire OSO (ou NULL si non disponible)
#' @param ndvi_mask SpatRaster binaire NDVI (ou NULL si non disponible)
#' @param method "union" ou "intersection"
#' @return SpatRaster binaire combiné (1 = forêt, 0 = non-forêt, NA = hors zone)
#' @export
build_forest_mask <- function(oso_mask = NULL, ndvi_mask = NULL,
                               method = FOREST_MASK_PARAMS$combine_method) {
  if (is.null(oso_mask) && is.null(ndvi_mask)) {
    return(NULL)
  }

  if (is.null(oso_mask)) return(ndvi_mask)
  if (is.null(ndvi_mask)) return(oso_mask)

  # Extraire les valeurs
  oso_vals  <- terra::values(oso_mask)[, 1]
  ndvi_vals <- terra::values(ndvi_mask)[, 1]

  combined <- if (method == "intersection") {
    as.integer(oso_vals == 1L & ndvi_vals == 1L)
  } else {
    # "union" par défaut
    as.integer(oso_vals == 1L | ndvi_vals == 1L)
  }

  r_combined <- terra::rast(oso_mask)
  terra::values(r_combined) <- combined
  names(r_combined) <- "forest_mask"

  r_combined
}

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

        # Créer un raster template aligné sur l'AOI à la résolution cible
        target_template <- terra::rast(aoi_ext, resolution = resolution, crs = "EPSG:2154")

        # Reprojeter directement vers le template (CRS + résolution + extent en une seule passe)
        r <- terra::project(r, target_template, method = "bilinear")

        # Masquer les pixels hors de l'AOI
        r <- terra::mask(r, aoi_vect)

        names(r) <- paste0(band, "_", format(current_date, "%Y%m%d"))
        band_stack[[band]] <- r
      }
    }

    if (length(band_stack) == length(bands)) {
      stacked <- terra::rast(band_stack)
      # terra::rast() peut perdre les noms custom → les remettre
      names(stacked) <- paste0(bands, "_", format(current_date, "%Y%m%d"))
      cube_list[[as.character(current_date)]] <- stacked
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
#' @param terrain_rasters Liste optionnelle de chemins raster terrain
#'   (sortie de compute_terrain_rasters : dem, slope, aspect, twi)
#' @param forest_mask SpatRaster binaire optionnel (1 = forêt, 0 = non-forêt).
#'   Si fourni, seuls les pixels forestiers sont traités (gain de temps + précision).
#' @return Liste avec la matrice de features, indices valides, etc.
extract_pixel_features <- function(cube_list, dates, block_size = 100,
                                    terrain_rasters = NULL,
                                    forest_mask = NULL) {
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
      layer_name <- grep(paste0("^", band, "(_|$)"), names(cube_list[[d]]), value = TRUE)
      if (length(layer_name) > 0) {
        vals <- terra::values(cube_list[[d]][[layer_name[1]]])
        band_mat[, d] <- as.numeric(vals)
      }
    }
    cube_arrays[[band]] <- band_mat
  }

  log_msg("  Cube chargé en mémoire")

  # Normalisation Sentinel-2 L2A : réflectance entière (0-10000) → 0-1
  # Les données synthétiques d'entraînement sont en réflectance 0-1,
  # les données S2 L2A réelles sont en réflectance × S2_SCALE_FACTOR
  first_max <- max(cube_arrays[[bands[1]]], na.rm = TRUE)
  if (!is.na(first_max) && first_max > 2) {
    log_msg("  Normalisation S2 : valeurs max={round(first_max)} → division par {S2_SCALE_FACTOR}")
    for (band in bands) {
      cube_arrays[[band]] <- cube_arrays[[band]] / S2_SCALE_FACTOR
    }
  }

  # Diagnostic : statistiques du cube pour la première bande
  first_band_mat <- cube_arrays[[bands[1]]]
  n_na_per_date <- colSums(is.na(first_band_mat))
  n_nonzero_per_date <- colSums(!is.na(first_band_mat) & first_band_mat != 0, na.rm = TRUE)
  log_msg("  Diagnostic {bands[1]} : {n_pixels} pixels × {n_dates} dates")
  log_msg("  NA par date (min/max) : {min(n_na_per_date)}/{max(n_na_per_date)}")
  log_msg("  Non-NA & non-zero par date (min/max) : {min(n_nonzero_per_date)}/{max(n_nonzero_per_date)}")

  if (all(is.na(first_band_mat))) {
    # Diagnostic supplémentaire : vérifier les noms de couches
    log_msg("  ATTENTION : toutes les valeurs sont NA !", level = "danger")
    log_msg("  Noms des couches du cube[1] : {paste(names(cube_list[[1]]), collapse=', ')}")
    log_msg("  Bandes recherchées : {paste(bands, collapse=', ')}")

    # Essayer de lire les valeurs brutes sans filtrage par nom
    test_vals <- terra::values(cube_list[[1]])
    n_not_na <- sum(!is.na(test_vals))
    log_msg("  Valeurs brutes non-NA dans cube[1] : {n_not_na} / {length(test_vals)}")
    if (n_not_na > 0) {
      val_range <- range(test_vals, na.rm = TRUE)
      log_msg("  Plage de valeurs : [{val_range[1]}, {val_range[2]}]")
    }
  }

  # --- Charger les rasters terrain si fournis ---
  terrain_vals <- NULL
  if (!is.null(terrain_rasters)) {
    log_msg("  Chargement des rasters terrain...")
    tryCatch({
      r_dem    <- terra::rast(terrain_rasters$dem)
      r_slope  <- terra::rast(terrain_rasters$slope)
      r_aspect <- terra::rast(terrain_rasters$aspect)
      r_twi    <- terra::rast(terrain_rasters$twi)
      r_tpi    <- terra::rast(terrain_rasters$tpi)

      # R\u00e9\u00e9chantillonner les rasters terrain sur la grille du cube S2
      r_dem    <- terra::resample(r_dem, template, method = "bilinear")
      r_slope  <- terra::resample(r_slope, template, method = "bilinear")
      r_aspect <- terra::resample(r_aspect, template, method = "bilinear")
      r_twi    <- terra::resample(r_twi, template, method = "bilinear")
      r_tpi    <- terra::resample(r_tpi, template, method = "bilinear")

      # Extraire toutes les valeurs en une passe (n_pixels \u00d7 1)
      dem_v    <- as.numeric(terra::values(r_dem))
      slope_v  <- as.numeric(terra::values(r_slope))
      aspect_v <- as.numeric(terra::values(r_aspect))
      twi_v    <- as.numeric(terra::values(r_twi))
      tpi_v    <- as.numeric(terra::values(r_tpi))

      # Convertir exposition en sin/cos
      aspect_rad <- aspect_v * pi / 180
      terrain_vals <- data.frame(
        DEM_elevation  = dem_v,
        DEM_slope      = slope_v,
        DEM_aspect_sin = sin(aspect_rad),
        DEM_aspect_cos = cos(aspect_rad),
        DEM_TWI        = twi_v,
        DEM_TPI        = tpi_v
      )
      log_msg("  Terrain charg\u00e9 : {nrow(terrain_vals)} pixels \u00d7 6 features", level = "success")
    }, error = function(e) {
      log_msg("  Erreur chargement terrain : {e$message}", level = "warning")
      terrain_vals <<- NULL
    })
  }

  # Identifier les pixels valides (au moins 3 dates non-NA pour la première bande)
  valid_mask <- rowSums(!is.na(first_band_mat)) >= 3
  n_data_valid <- sum(valid_mask)

  # --- Appliquer le masque forestier (OSO + NDVI) ---
  if (!is.null(forest_mask)) {
    forest_vals <- as.integer(terra::values(forest_mask)[, 1])
    forest_bool <- !is.na(forest_vals) & forest_vals == 1L
    n_forest <- sum(forest_bool)
    n_masked_out <- sum(valid_mask & !forest_bool)
    valid_mask <- valid_mask & forest_bool
    log_msg("  Masque forestier : {n_forest} pixels forêt, {n_masked_out} pixels non-forêt exclus")
  }

  valid_idx  <- which(valid_mask)
  n_valid    <- length(valid_idx)

  log_msg("  {n_valid} pixels valides sur {n_pixels} ({round(n_valid/n_pixels*100,1)}%)")

  if (n_valid == 0) {
    cli::cli_alert_danger("Aucun pixel valide dans l'AOI")
    cli::cli_alert_info("Vérifiez que les rasters téléchargés couvrent bien la zone d'intérêt.")
    return(NULL)
  }

  # --- Extraction des features pour chaque pixel valide ---
  log_msg("  Extraction des features...")

  # Labels DOY pour les dates cibles (cohérent avec build_feature_matrix)
  doy_labels <- format(target_dates, "%j")

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

    # 1. Bandes brutes : série temporelle + statistiques
    for (band in bands) {
      vals <- bands_ts[[band]]

      # Série temporelle brute par DOY (comme dans build_feature_matrix)
      feat_names <- paste0(band, "_d", doy_labels)
      features <- c(features, setNames(vals, feat_names))

      # Statistiques temporelles
      stats <- calc_temporal_stats(vals, prefix = band)
      features <- c(features, stats)
    }

    # 2. Indices spectraux : série temporelle + statistiques
    if (all(c("B08", "B04", "B02", "B03", "B05", "B11", "B12") %in% names(bands_ts))) {
      indices <- calc_all_indices(bands_ts)

      for (idx_name in names(indices)) {
        idx_vals <- indices[[idx_name]]

        # Série temporelle brute par DOY
        idx_feat_names <- paste0(idx_name, "_d", doy_labels[seq_along(idx_vals)])
        features <- c(features, setNames(idx_vals, idx_feat_names))

        # Statistiques temporelles
        stats <- calc_temporal_stats(idx_vals, prefix = idx_name)
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

    # 5. Features terrain (MNT, pente, exposition sin/cos, TWI, TPI)
    if (!is.null(terrain_vals)) {
      px_terrain <- c(
        DEM_elevation  = terrain_vals$DEM_elevation[px_idx],
        DEM_slope      = terrain_vals$DEM_slope[px_idx],
        DEM_aspect_sin = terrain_vals$DEM_aspect_sin[px_idx],
        DEM_aspect_cos = terrain_vals$DEM_aspect_cos[px_idx],
        DEM_TWI        = terrain_vals$DEM_TWI[px_idx],
        DEM_TPI        = terrain_vals$DEM_TPI[px_idx]
      )
      # Remplacer les NA terrain par 0
      px_terrain[is.na(px_terrain)] <- 0
      features <- c(features, px_terrain)
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

#' Reconstruction des rasters classifiés à partir des prédictions pixellaires
#'
#' Produit une classification soft complète :
#' - Raster catégoriel de l'essence dominante (comme avant)
#' - Raster de confiance (probabilité max)
#' - Raster multi-bandes des probabilités par essence (1 bande = 1 essence)
#' - Raster d'entropie de Shannon (taux de mélange)
#' - Raster multi-bandes de présence binaire par essence (proba > seuil)
#'
#' @param predictions Sortie de classify_pixels (contient all_probas)
#' @param pixel_features Sortie de extract_pixel_features
#' @return Liste avec tous les rasters et métadonnées
build_species_raster <- function(predictions, pixel_features) {
  log_msg("Construction des rasters (classification soft)...")

  template <- pixel_features$template
  n_pixels <- pixel_features$n_rows * pixel_features$n_cols
  class_names <- predictions$class_names
  n_classes <- length(class_names)
  valid_idx <- predictions$valid_idx

  # --- 1. Raster des classes (essence dominante, code 1-N) ---
  class_vals <- rep(NA_real_, n_pixels)
  class_vals[valid_idx] <- predictions$class_idx

  r_class <- terra::rast(template)
  terra::values(r_class) <- class_vals
  names(r_class) <- "species_code"

  levels_df <- data.frame(
    value   = seq_along(class_names),
    species = class_names
  )
  levels(r_class) <- levels_df

  # --- 2. Raster de confiance (probabilité max) ---
  proba_vals <- rep(NA_real_, n_pixels)
  proba_vals[valid_idx] <- predictions$max_proba

  r_confidence <- terra::rast(template)
  terra::values(r_confidence) <- proba_vals
  names(r_confidence) <- "confidence"

  # --- 3. Raster multi-bandes : probabilités par essence ---
  r_probas <- NULL
  if (isTRUE(CLASSIF_PARAMS$export_probabilities)) {
    log_msg("  Construction raster multi-bandes probabilités ({n_classes} bandes)...")
    proba_layers <- list()
    for (k in seq_len(n_classes)) {
      vals <- rep(NA_real_, n_pixels)
      vals[valid_idx] <- predictions$all_probas[, k]
      r_k <- terra::rast(template)
      terra::values(r_k) <- vals
      proba_layers[[k]] <- r_k
    }
    r_probas <- terra::rast(proba_layers)
    names(r_probas) <- paste0("prob_", gsub(" ", "_", class_names))
  }

  # --- 4. Raster entropie de Shannon (taux de mélange) ---
  r_shannon <- NULL
  if (isTRUE(CLASSIF_PARAMS$export_shannon)) {
    log_msg("  Calcul de l'entropie de Shannon par pixel...")
    # H = -sum(p * log(p)), normalisée par log(n_classes) → [0, 1]
    shannon_vals_valid <- apply(predictions$all_probas, 1, function(p) {
      p <- p[p > 0]  # éviter log(0)
      -sum(p * log(p)) / log(n_classes)
    })

    shannon_vals <- rep(NA_real_, n_pixels)
    shannon_vals[valid_idx] <- shannon_vals_valid

    r_shannon <- terra::rast(template)
    terra::values(r_shannon) <- shannon_vals
    names(r_shannon) <- "shannon_entropy"
  }

  # --- 5. Raster multi-bandes : présence binaire (proba > seuil) ---
  r_presence <- NULL
  if (isTRUE(CLASSIF_PARAMS$export_presence)) {
    threshold <- CLASSIF_PARAMS$presence_threshold %||% 0.10
    log_msg("  Cartes de présence par essence (seuil = {threshold})...")
    presence_layers <- list()
    for (k in seq_len(n_classes)) {
      vals <- rep(NA_integer_, n_pixels)
      vals[valid_idx] <- as.integer(predictions$all_probas[, k] >= threshold)
      r_k <- terra::rast(template)
      terra::values(r_k) <- vals
      presence_layers[[k]] <- r_k
    }
    r_presence <- terra::rast(presence_layers)
    names(r_presence) <- paste0("presence_", gsub(" ", "_", class_names))
  }

  log_msg("  Rasters construits (soft classification)", level = "success")

  list(
    species    = r_class,
    confidence = r_confidence,
    probas     = r_probas,
    shannon    = r_shannon,
    presence   = r_presence,
    legend     = levels_df
  )
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
                                 use_s1       = TRUE,
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
        models_dir <- file.path(.get_project_root(), "output", "models")
        default_pt <- file.path(models_dir, paste0("treesatai_", pytorch_model, "_best.pt"))
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
    models_dir <- file.path(.get_project_root(), "output", "models")
    default_rds <- file.path(models_dir, "treesatai_rf.rds")

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

  # --- 3b. Données terrain (MNT, pente, exposition, TWI, TPI) ---
  terrain_rasters <- NULL
  if (isTRUE(CLASSIF_PARAMS$use_terrain)) {
    cli::cli_h2("3b. Données terrain (MNT + dérivés)")

    dem_dir <- file.path(RAW_DIR, "dem")
    dir.create(dem_dir, showWarnings = FALSE, recursive = TRUE)

    # Déterminer la source DEM
    # - "auto" : détecte si l'AOI est en France → IGN 1 m, sinon Copernicus 30 m
    # - "ign"  : force IGN (France uniquement)
    # - "copernicus" : force Copernicus (Europe entière)
    dem_src <- DEM_PARAMS$dem_source
    if (dem_src == "auto") {
      # Vérifier si l'AOI est en France métropolitaine (bbox approx)
      aoi_wgs <- sf::st_transform(aoi, 4326)
      bbox_wgs <- sf::st_bbox(aoi_wgs)
      in_france <- bbox_wgs["xmin"] >= -5.5 && bbox_wgs["xmax"] <= 10 &&
                   bbox_wgs["ymin"] >= 41 && bbox_wgs["ymax"] <= 51.5
      dem_src <- if (in_france) "ign" else "copernicus"
      log_msg("  Source DEM auto-détectée : {dem_src}")
    }

    # Télécharger le MNT
    dem_path <- tryCatch({
      if (dem_src == "ign") {
        download_dem_ign(aoi, output_dir = dem_dir,
                         resolution = DEM_PARAMS$resample_res)
      } else {
        download_dem_copernicus(aoi, output_dir = dem_dir)
      }
    }, error = function(e) {
      log_msg("Erreur téléchargement MNT ({dem_src}) : {e$message}", level = "warning")
      # Fallback : si IGN échoue, essayer Copernicus
      if (dem_src == "ign") {
        log_msg("Tentative fallback Copernicus DEM 30 m...", level = "warning")
        tryCatch(
          download_dem_copernicus(aoi, output_dir = dem_dir),
          error = function(e2) {
            log_msg("Impossible d'obtenir un MNT : {e2$message}", level = "danger")
            NULL
          }
        )
      } else {
        NULL
      }
    })

    if (!is.null(dem_path) && file.exists(dem_path)) {
      terrain_rasters <- tryCatch(
        compute_terrain_rasters(dem_path, output_dir = dem_dir),
        error = function(e) {
          log_msg("Erreur calcul dérivés terrain : {e$message}", level = "warning")
          NULL
        }
      )
      if (!is.null(terrain_rasters)) {
        log_msg("Terrain prêt : MNT + pente + exposition + TWI + TPI", level = "success")
      }
    } else {
      log_msg("Pas de MNT disponible — classification sans features terrain", level = "warning")
    }
  }

  # --- 3c. Masque forestier (OSO + NDVI) ---
  forest_mask_raster <- NULL
  if (isTRUE(FOREST_MASK_PARAMS$apply_forest_mask)) {
    cli::cli_h2("3c. Masque forestier (OSO + NDVI)")

    template <- cube_list[[1]][[1]]
    oso_mask <- NULL
    ndvi_mask <- NULL

    # --- OSO ---
    if (isTRUE(FOREST_MASK_PARAMS$use_oso)) {
      oso_path <- tryCatch({
        download_oso(aoi, year = FOREST_MASK_PARAMS$oso_year)
      }, error = function(e) {
        log_msg("Erreur OSO : {e$message}", level = "warning")
        NULL
      })

      if (!is.null(oso_path) && file.exists(oso_path)) {
        oso_mask <- build_oso_forest_mask(oso_path, template)
        n_oso_forest <- sum(terra::values(oso_mask) == 1L, na.rm = TRUE)
        n_total <- terra::ncell(oso_mask)
        log_msg("  OSO : {n_oso_forest}/{n_total} pixels forestiers ({round(n_oso_forest/n_total*100,1)}%)",
                level = "success")
      } else {
        log_msg("  OSO non disponible — masque NDVI seul", level = "warning")
      }
    }

    # --- NDVI ---
    if (isTRUE(FOREST_MASK_PARAMS$use_ndvi)) {
      log_msg("  Calcul du masque NDVI (seuil = {FOREST_MASK_PARAMS$ndvi_min_threshold})...")

      # Charger les bandes B08 et B04 du cube pour le calcul NDVI
      bands <- S2_BAND_NAMES
      n_dates_cube <- length(cube_list)
      n_pixels_cube <- terra::ncell(template)
      cube_b08 <- matrix(NA_real_, nrow = n_pixels_cube, ncol = n_dates_cube)
      cube_b04 <- matrix(NA_real_, nrow = n_pixels_cube, ncol = n_dates_cube)
      for (d in seq_len(n_dates_cube)) {
        b08_name <- grep("^B08(_|$)", names(cube_list[[d]]), value = TRUE)
        b04_name <- grep("^B04(_|$)", names(cube_list[[d]]), value = TRUE)
        if (length(b08_name) > 0) cube_b08[, d] <- as.numeric(terra::values(cube_list[[d]][[b08_name[1]]]))
        if (length(b04_name) > 0) cube_b04[, d] <- as.numeric(terra::values(cube_list[[d]][[b04_name[1]]]))
      }

      # Normaliser si nécessaire (valeurs S2 L2A 0-10000 → 0-1)
      if (max(cube_b08, na.rm = TRUE) > 2) {
        cube_b08 <- cube_b08 / S2_SCALE_FACTOR
        cube_b04 <- cube_b04 / S2_SCALE_FACTOR
      }

      ndvi_arrays <- list(B08 = cube_b08, B04 = cube_b04)
      ndvi_mask <- build_ndvi_mask(ndvi_arrays, template)
      n_ndvi_forest <- sum(terra::values(ndvi_mask) == 1L, na.rm = TRUE)
      log_msg("  NDVI : {n_ndvi_forest}/{n_pixels_cube} pixels végétation active",
              level = "success")
    }

    # --- Combiner ---
    forest_mask_raster <- build_forest_mask(oso_mask, ndvi_mask)
    if (!is.null(forest_mask_raster)) {
      n_forest <- sum(terra::values(forest_mask_raster) == 1L, na.rm = TRUE)
      n_total <- terra::ncell(forest_mask_raster)
      log_msg("  Masque combiné ({FOREST_MASK_PARAMS$combine_method}) : {n_forest}/{n_total} pixels forestiers ({round(n_forest/n_total*100,1)}%)",
              level = "success")

      # Sauvegarder le masque pour inspection
      mask_path <- file.path(output_dir, "masque_foret.tif")
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      if (file.exists(mask_path)) file.remove(mask_path)
      terra::writeRaster(forest_mask_raster, mask_path, datatype = "INT1U")
      log_msg("  Masque sauvegardé : {mask_path}", level = "info")
    }
  }

  # --- 4. Extraction des features ---
  cli::cli_h2("4. Extraction des features pixellaires")
  pixel_features <- extract_pixel_features(cube_list, dates,
                                            terrain_rasters = terrain_rasters,
                                            forest_mask = forest_mask_raster)

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
  output_files <- list()

  # 7a. Raster essence dominante (catégoriel, comme avant)
  tif_path <- file.path(output_dir, "carte_essences.tif")
  if (file.exists(tif_path)) file.remove(tif_path)
  terra::writeRaster(rasters$species, tif_path, datatype = "INT1U")
  log_msg("  Raster espèces       : {tif_path}", level = "success")
  output_files$species <- tif_path

  # 7b. Raster confiance (probabilité max)
  conf_path <- file.path(output_dir, "carte_confiance.tif")
  if (file.exists(conf_path)) file.remove(conf_path)
  terra::writeRaster(rasters$confidence, conf_path)
  log_msg("  Raster confiance     : {conf_path}", level = "success")
  output_files$confidence <- conf_path

  # 7c. Raster multi-bandes probabilités (1 bande par essence)
  if (!is.null(rasters$probas)) {
    probas_path <- file.path(output_dir, "carte_probabilites.tif")
    if (file.exists(probas_path)) file.remove(probas_path)
    terra::writeRaster(rasters$probas, probas_path)
    log_msg("  Raster probabilités  : {probas_path} ({terra::nlyr(rasters$probas)} bandes)",
            level = "success")
    output_files$probas <- probas_path
  }

  # 7d. Raster entropie de Shannon (taux de mélange, 0 = pur, 1 = mélange max)
  if (!is.null(rasters$shannon)) {
    shannon_path <- file.path(output_dir, "carte_shannon.tif")
    if (file.exists(shannon_path)) file.remove(shannon_path)
    terra::writeRaster(rasters$shannon, shannon_path)
    log_msg("  Raster Shannon       : {shannon_path}", level = "success")
    output_files$shannon <- shannon_path
  }

  # 7e. Raster multi-bandes présence binaire (proba >= seuil)
  if (!is.null(rasters$presence)) {
    presence_path <- file.path(output_dir, "carte_presence.tif")
    if (file.exists(presence_path)) file.remove(presence_path)
    terra::writeRaster(rasters$presence, presence_path, datatype = "INT1U")
    log_msg("  Raster présence      : {presence_path} ({terra::nlyr(rasters$presence)} bandes)",
            level = "success")
    output_files$presence <- presence_path
  }

  # 7f. Vectoriser le raster essence dominante (polygones)
  log_msg("  Vectorisation...")
  species_poly <- terra::as.polygons(rasters$species, dissolve = TRUE)
  species_sf   <- sf::st_as_sf(species_poly)

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
  log_msg("  Vecteur espèces      : {gpkg_path}", level = "success")
  output_files$vector <- gpkg_path

  # 7g. Légende CSV
  legend_path <- file.path(output_dir, "legende_especes.csv")
  readr::write_csv(rasters$legend, legend_path)

  # 7h. Statistiques
  stats <- compute_map_statistics(predictions, rasters)
  stats_path <- file.path(output_dir, "statistiques_essences.csv")
  readr::write_csv(stats, stats_path)
  log_msg("  Statistiques         : {stats_path}", level = "success")

  # 7i. Rapport cartographique PDF
  s2_rgb <- tryCatch(
    build_s2_rgb_composite(cube_list, dates),
    error = function(e) {
      log_msg("  Composite RGB : {e$message}", level = "warning")
      NULL
    }
  )
  # 7j. Rapport cartographique (dashboard mono-page : RStudio + PDF)
  report <- tryCatch(
    generate_prediction_report(
      rasters     = rasters,
      statistics  = stats,
      output_dir  = output_dir,
      s2_rgb      = s2_rgb,
      aoi         = aoi,
      forest_mask = forest_mask_raster
    ),
    error = function(e) {
      log_msg("  Erreur rapport cartographique : {e$message}", level = "warning")
      NULL
    }
  )
  if (!is.null(report)) {
    output_files$dashboard  <- report$dashboard
    output_files$pdf_report <- report$pdf_path
  }

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
  cli::cli_li("{output_files$species}    — essence dominante (GeoTIFF catégoriel)")
  cli::cli_li("{output_files$confidence} — confiance (probabilité max)")
  if (!is.null(output_files$probas))
    cli::cli_li("{output_files$probas}   — probabilités par essence (multi-bandes)")
  if (!is.null(output_files$shannon))
    cli::cli_li("{output_files$shannon}  — entropie de Shannon (taux de mélange)")
  if (!is.null(output_files$presence))
    cli::cli_li("{output_files$presence} — présence par essence (binaire, multi-bandes)")
  cli::cli_li("{output_files$vector}     — polygones par espèce (GeoPackage)")
  cli::cli_li("{stats_path}              — statistiques par espèce")
  if (!is.null(forest_mask_raster))
    cli::cli_li("{file.path(output_dir, 'masque_foret.tif')} — masque forestier (OSO+NDVI)")
  if (!is.null(output_files$pdf_report))
    cli::cli_li("{output_files$pdf_report} — rapport cartographique (PDF 8 pages)")
  if (!is.null(output_files$pdf_rstudio))
    cli::cli_li("{output_files$pdf_rstudio} — rapport RStudio (ggplot2+patchwork)")
  cli::cli_end()

  if (!is.null(output_files$shannon)) {
    cli::cli_text("")
    cli::cli_alert_info(paste0(
      "Entropie de Shannon normalisée [0-1] : ",
      "0 = peuplement pur (100% une essence), ",
      "1 = mélange maximal (équiprobable). ",
      "Utilisable directement pour l'indice de biodiversité Néméton (famille B)."
    ))
  }

  invisible(list(
    species_raster    = rasters$species,
    confidence_raster = rasters$confidence,
    probas_raster     = rasters$probas,
    shannon_raster    = rasters$shannon,
    presence_raster   = rasters$presence,
    forest_mask       = forest_mask_raster,
    species_vector    = species_sf,
    statistics        = stats,
    model             = model,
    legend            = rasters$legend,
    pdf_report        = pdf_report,
    rstudio_report    = rstudio_report
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

  # Dates cibles pour les labels DOY (cohérent avec build_feature_matrix)
  s1_dates <- as.Date(names(s1_cube))
  s1_year <- format(s1_dates[1], "%Y")
  s1_target_dates <- seq.Date(
    as.Date(paste0(s1_year, "-01-01")),
    as.Date(paste0(s1_year, "-12-31")),
    by = TS_PARAMS$target_interval_days
  )
  s1_doy_labels <- format(s1_target_dates, "%j")

  # Extraire les features S1 par pixel
  feature_list <- vector("list", n_valid)

  for (i in seq_len(n_valid)) {
    px <- valid_idx[i]
    vv_raw <- s1_arrays$VV[px, ]
    vh_raw <- s1_arrays$VH[px, ]

    # Interpolation aux dates cibles (comme pour S2)
    vv_ts <- interpolate_ts(s1_dates, vv_raw, s1_target_dates)
    vh_ts <- interpolate_ts(s1_dates, vh_raw, s1_target_dates)

    # Lissage
    if (!all(is.na(vv_ts))) vv_ts <- smooth_savgol(vv_ts)
    if (!all(is.na(vh_ts))) vh_ts <- smooth_savgol(vh_ts)

    features <- c()

    # Série temporelle brute par DOY (cohérent avec build_feature_matrix)
    features <- c(features, setNames(vv_ts, paste0("S1_VV_d", s1_doy_labels)))
    features <- c(features, setNames(vh_ts, paste0("S1_VH_d", s1_doy_labels)))

    # Statistiques temporelles + indices radar
    features <- c(features, calc_s1_temporal_features(vv_ts, vh_ts))

    feature_list[[i]] <- features
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

  # Entropie de Shannon moyenne par pixel dominant de chaque espèce
  if (!is.null(rasters$shannon)) {
    shannon_valid <- terra::values(rasters$shannon)[predictions$valid_idx]
    shannon_by_class <- tapply(shannon_valid, class_idx, mean, na.rm = TRUE)
    stats$shannon_moy <- round(as.numeric(shannon_by_class[as.character(stats$code)]), 3)
    stats$shannon_moy[is.na(stats$shannon_moy)] <- 0
  }

  # Trier par surface décroissante
  stats <- stats[order(-stats$n_pixels), ]

  # Ajouter le type (compatible avec les 10 groupes ou les 20 espèces)
  if (all(stats$espece %in% SPECIES_GROUPS_INFO$group)) {
    stats <- stats |>
      dplyr::left_join(SPECIES_GROUPS_INFO[, c("group", "type", "phenologie")],
                       by = c("espece" = "group"))
  } else {
    stats <- stats |>
      dplyr::left_join(SPECIES[, c("french", "type", "phenologie")],
                       by = c("espece" = "french"))
  }

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

