#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Fonctions utilitaires
# ==============================================================================

# --- Calcul des indices spectraux Sentinel-2 ----------------------------------

#' Calcul du NDVI (Normalized Difference Vegetation Index)
#' @param nir Bande NIR (B08)
#' @param red Bande Rouge (B04)
calc_ndvi <- function(nir, red) {
  (nir - red) / (nir + red + 1e-10)
}

#' Calcul de l'EVI (Enhanced Vegetation Index)
#' @param nir Bande NIR (B08)
#' @param red Bande Rouge (B04)
#' @param blue Bande Bleue (B02)
calc_evi <- function(nir, red, blue) {
  2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1 + 1e-10)
}

#' Calcul du NDWI (Normalized Difference Water Index)
#' @param nir Bande NIR (B08)
#' @param swir Bande SWIR1 (B11)
calc_ndwi <- function(nir, swir) {
  (nir - swir) / (nir + swir + 1e-10)
}

#' Calcul du CRI (Carotenoid Reflectance Index)
#' @param green Bande verte (B03)
#' @param rededge1 Bande Red Edge 1 (B05)
calc_cri <- function(green, rededge1) {
  (1 / (green + 1e-10)) - (1 / (rededge1 + 1e-10))
}

#' Calcul du RENDVI (Red Edge NDVI)
#' @param nir Bande NIR (B08)
#' @param rededge1 Bande Red Edge 1 (B05)
calc_rendvi <- function(nir, rededge1) {
  (nir - rededge1) / (nir + rededge1 + 1e-10)
}

#' Calcul du NBR (Normalized Burn Ratio)
#' @param nir Bande NIR (B08)
#' @param swir2 Bande SWIR2 (B12)
calc_nbr <- function(nir, swir2) {
  (nir - swir2) / (nir + swir2 + 1e-10)
}

#' Calcul de tous les indices spectraux pour un pixel/parcelle
#' @param bands Liste nommée des bandes spectrales (vecteurs temporels)
#' @return data.frame avec les indices spectraux par date
calc_all_indices <- function(bands) {
  data.frame(
    NDVI   = calc_ndvi(bands$B08, bands$B04),
    EVI    = calc_evi(bands$B08, bands$B04, bands$B02),
    NDWI   = calc_ndwi(bands$B08, bands$B11),
    CRI    = calc_cri(bands$B03, bands$B05),
    RENDVI = calc_rendvi(bands$B08, bands$B05),
    NBR    = calc_nbr(bands$B08, bands$B12)
  )
}

# --- Masquage nuageux --------------------------------------------------------

#' Application du masque nuageux SCL (Scene Classification Layer) Sentinel-2
#' @param scl Valeurs SCL
#' @return Vecteur logique TRUE = pixel valide (clair)
scl_clear_mask <- function(scl) {
  # SCL classes valides :
  #  4 = Vegetation
  #  5 = Bare soils
  #  6 = Water
  #  7 = Unclassified (on le garde prudemment)
  # Classes masquées :
  #  0 = No data, 1 = Saturated, 2 = Dark area, 3 = Cloud shadow,
  #  8 = Cloud medium, 9 = Cloud high, 10 = Thin cirrus, 11 = Snow
  scl %in% c(4, 5, 6)
}

# --- Interpolation temporelle ------------------------------------------------

#' Interpolation linéaire des séries temporelles avec données manquantes
#' @param dates Vecteur de dates d'acquisition
#' @param values Vecteur de valeurs (avec NA pour nuages)
#' @param target_dates Vecteur de dates cibles pour l'interpolation
#' @return Vecteur interpolé aux dates cibles
interpolate_ts <- function(dates, values, target_dates) {
  # Retirer les NA
  valid <- !is.na(values)
  if (sum(valid) < 3) {
    return(rep(NA_real_, length(target_dates)))
  }

  # Convertir en numérique (jours depuis le début)
  origin <- min(target_dates)
  x_valid  <- as.numeric(difftime(dates[valid], origin, units = "days"))
  y_valid  <- values[valid]
  x_target <- as.numeric(difftime(target_dates, origin, units = "days"))

  # Interpolation par spline
  tryCatch({
    spline_fit <- stats::spline(x_valid, y_valid, xout = x_target, method = "natural")
    spline_fit$y
  }, error = function(e) {
    # Fallback sur interpolation linéaire
    stats::approx(x_valid, y_valid, xout = x_target, rule = 2)$y
  })
}

# --- Lissage Savitzky-Golay ---------------------------------------------------

#' Lissage d'une série temporelle par filtre Savitzky-Golay
#' @param ts_values Vecteur de valeurs temporelles
#' @param order Ordre du polynôme (défaut 3)
#' @param length Taille de la fenêtre (impair, défaut 7)
#' @return Vecteur lissé
smooth_savgol <- function(ts_values, order = 3, length = 7) {
  if (any(is.na(ts_values))) {
    # Remplir les NA par interpolation linéaire d'abord
    idx <- seq_along(ts_values)
    valid <- !is.na(ts_values)
    if (sum(valid) < length) return(ts_values)
    ts_values <- stats::approx(idx[valid], ts_values[valid], xout = idx, rule = 2)$y
  }

  if (requireNamespace("signal", quietly = TRUE)) {
    # Filtre Savitzky-Golay du package signal
    signal::sgolayfilt(ts_values, p = order, n = length)
  } else {
    # Fallback : moyenne mobile pondérée
    stats::filter(ts_values, rep(1 / length, length), sides = 2) |>
      as.numeric() |>
      (\(x) ifelse(is.na(x), ts_values, x))()
  }
}

# --- Normalisation ------------------------------------------------------------

#' Normalisation Min-Max d'un vecteur
#' @param x Vecteur numérique
#' @return Vecteur normalisé entre 0 et 1
normalize_minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] == rng[1]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

#' Normalisation Z-score
#' @param x Vecteur numérique
#' @return Vecteur centré-réduit
normalize_zscore <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

# --- Statistiques temporelles -------------------------------------------------

#' Calcul de statistiques descriptives sur une série temporelle
#' @param ts_values Vecteur de valeurs temporelles
#' @param prefix Préfixe pour les noms des features
#' @return Vecteur nommé de statistiques
calc_temporal_stats <- function(ts_values, prefix = "NDVI") {
  ts_clean <- ts_values[!is.na(ts_values)]
  if (length(ts_clean) == 0) {
    return(setNames(rep(NA_real_, 9), paste0(prefix, "_",
      c("mean", "sd", "min", "max", "range", "median", "q25", "q75", "iqr"))))
  }

  q <- quantile(ts_clean, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

  stats <- c(
    mean   = mean(ts_clean, na.rm = TRUE),
    sd     = sd(ts_clean, na.rm = TRUE),
    min    = min(ts_clean, na.rm = TRUE),
    max    = max(ts_clean, na.rm = TRUE),
    range  = diff(range(ts_clean, na.rm = TRUE)),
    median = q[2],
    q25    = q[1],
    q75    = q[3],
    iqr    = q[3] - q[1]
  )
  setNames(stats, paste0(prefix, "_", names(stats)))
}

# --- Gestion des fichiers ----------------------------------------------------

#' Liste les fichiers raster Sentinel-2 par bande et date
#' @param input_dir Répertoire contenant les images Sentinel-2
#' @param band Nom de la bande (ex: "B04")
#' @return data.frame avec colonnes : path, date, band
list_s2_files <- function(input_dir, band = NULL) {
  pattern <- if (!is.null(band)) {
    glue::glue(".*{band}.*\\.(tif|TIF|jp2|JP2)$")
  } else {
    ".*\\.(tif|TIF|jp2|JP2)$"
  }

  files <- list.files(input_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)

  if (length(files) == 0) {
    cli::cli_warn("Aucun fichier trouvé dans {input_dir} pour la bande {band}")
    return(data.frame(path = character(), date = as.Date(character()), band = character()))
  }

  # Extraction de la date depuis le nom de fichier Sentinel-2
  # Format typique : S2A_MSIL2A_20210315T105031_...
  date_pattern <- "(\\d{4})(\\d{2})(\\d{2})"
  dates <- stringr::str_extract(basename(files), date_pattern)
  parsed_dates <- as.Date(dates, format = "%Y%m%d")

  data.frame(
    path = files,
    date = parsed_dates,
    band = if (!is.null(band)) band else stringr::str_extract(basename(files), "B\\d{2}|B8A"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(date)
}

# --- Logging ------------------------------------------------------------------

#' Log un message avec horodatage
log_msg <- function(..., level = "info") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(...)
  switch(level,
    "info"    = cli::cli_alert_info("[{timestamp}] {msg}"),
    "success" = cli::cli_alert_success("[{timestamp}] {msg}"),
    "warning" = cli::cli_alert_warning("[{timestamp}] {msg}"),
    "danger"  = cli::cli_alert_danger("[{timestamp}] {msg}"),
    cli::cli_text("[{timestamp}] {msg}")
  )
}
