#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Extraction de métriques phénologiques
# Signatures phénologiques annuelles pour 20 espèces européennes
# ==============================================================================

# Fonctions utilitaires chargées via le package (01_utils.R)

# --- Détection des seuils phénologiques ---------------------------------------

#' Détection du SOS (Start of Season) et EOS (End of Season) par seuillage
#' @param ndvi_ts Vecteur NDVI lissé (série temporelle annuelle)
#' @param dates Vecteur de dates correspondantes
#' @param threshold Seuil relatif (proportion de l'amplitude, défaut 0.5)
#' @return Liste avec SOS, EOS, et informations associées
detect_season <- function(ndvi_ts, dates, threshold = 0.5) {
  if (all(is.na(ndvi_ts))) {
    return(list(
      SOS = NA, EOS = NA, LOS = NA,
      SOS_value = NA, EOS_value = NA,
      max_date = NA, max_value = NA
    ))
  }

  doy <- as.numeric(format(dates, "%j"))

  # Valeurs min et max
  ndvi_min  <- min(ndvi_ts, na.rm = TRUE)
  ndvi_max  <- max(ndvi_ts, na.rm = TRUE)
  amplitude <- ndvi_max - ndvi_min

  # Seuil absolu
  thresh_value <- ndvi_min + threshold * amplitude

  # Date du maximum
  max_idx  <- which.max(ndvi_ts)
  max_date <- dates[max_idx]
  max_doy  <- doy[max_idx]

  # SOS : première date avant le max où NDVI dépasse le seuil (en montant)
  ascending <- ndvi_ts[1:max_idx]
  sos_idx <- which(ascending >= thresh_value)[1]
  SOS <- if (!is.na(sos_idx)) doy[sos_idx] else NA

  # EOS : première date après le max où NDVI repasse sous le seuil (en descendant)
  if (max_idx < length(ndvi_ts)) {
    descending <- ndvi_ts[max_idx:length(ndvi_ts)]
    eos_rel_idx <- which(descending <= thresh_value)[1]
    eos_idx <- if (!is.na(eos_rel_idx)) max_idx + eos_rel_idx - 1 else NA
    EOS <- if (!is.na(eos_idx)) doy[eos_idx] else NA
  } else {
    eos_idx <- NA
    EOS <- NA
  }

  # LOS : Length of Season
  LOS <- if (!is.na(SOS) && !is.na(EOS)) EOS - SOS else NA

  list(
    SOS       = SOS,
    EOS       = EOS,
    LOS       = LOS,
    SOS_value = if (!is.na(sos_idx)) ndvi_ts[sos_idx] else NA,
    EOS_value = if (!is.na(eos_idx)) ndvi_ts[eos_idx] else NA,
    max_date  = max_doy,
    max_value = ndvi_max
  )
}

# --- Extraction complète des métriques phénologiques -------------------------

#' Extraction d'un ensemble complet de métriques phénologiques
#' @param ndvi_ts Vecteur NDVI lissé
#' @param dates Vecteur de dates
#' @return Vecteur nommé de métriques phénologiques
extract_phenometrics <- function(ndvi_ts, dates) {
  if (all(is.na(ndvi_ts))) {
    return(setNames(rep(NA_real_, 28), paste0("pheno_", 1:28)))
  }

  doy <- as.numeric(format(dates, "%j"))
  n <- length(ndvi_ts)

  # Retirer les NA de façon synchrone (trapz exige des vecteurs sans NA)
  valid <- !is.na(ndvi_ts) & !is.na(doy)
  dates   <- dates[valid]
  doy     <- doy[valid]
  ndvi_ts <- ndvi_ts[valid]
  n <- length(ndvi_ts)

  if (n < 3) {
    return(setNames(rep(NA_real_, 28), paste0("pheno_", 1:28)))
  }

  # 1. Saisonnalité basique
  season <- detect_season(ndvi_ts, dates, threshold = 0.5)

  # 2. Amplitude saisonnière
  ndvi_min <- min(ndvi_ts, na.rm = TRUE)
  ndvi_max <- max(ndvi_ts, na.rm = TRUE)
  amplitude <- ndvi_max - ndvi_min

  # 3. Intégrales (proxy de productivité)
  # Intégrale totale sur l'année
  total_integral <- pracma::trapz(doy, ndvi_ts)

  # Intégrale de la saison de croissance (entre SOS et EOS)
  if (!is.na(season$SOS) && !is.na(season$EOS)) {
    growing_mask <- doy >= season$SOS & doy <= season$EOS
    if (sum(growing_mask) > 1) {
      growing_integral <- pracma::trapz(doy[growing_mask], ndvi_ts[growing_mask])
    } else {
      growing_integral <- NA_real_
    }
  } else {
    growing_integral <- NA_real_
  }

  # 4. Taux de verdissement (greening rate) — pente du NDVI au printemps
  spring_mask <- doy >= 60 & doy <= 180  # Mars à Juin
  if (sum(spring_mask) > 2) {
    spring_fit <- lm(ndvi_ts[spring_mask] ~ doy[spring_mask])
    greening_rate <- coef(spring_fit)[2]
  } else {
    greening_rate <- NA_real_
  }

  # 5. Taux de sénescence (browning rate) — pente du NDVI en automne
  autumn_mask <- doy >= 240 & doy <= 340  # Sept à Déc
  if (sum(autumn_mask) > 2) {
    autumn_fit <- lm(ndvi_ts[autumn_mask] ~ doy[autumn_mask])
    browning_rate <- coef(autumn_fit)[2]
  } else {
    browning_rate <- NA_real_
  }

  # 6. NDVI par saison (moyenne saisonnière)
  winter_mask <- doy <= 80 | doy >= 335    # Déc-Mars
  spring2_mask <- doy > 80 & doy <= 172    # Mars-Juin
  summer_mask <- doy > 172 & doy <= 264    # Juin-Sept
  autumn2_mask <- doy > 264 & doy < 335    # Sept-Déc

  ndvi_winter <- mean(ndvi_ts[winter_mask], na.rm = TRUE)
  ndvi_spring <- mean(ndvi_ts[spring2_mask], na.rm = TRUE)
  ndvi_summer <- mean(ndvi_ts[summer_mask], na.rm = TRUE)
  ndvi_autumn <- mean(ndvi_ts[autumn2_mask], na.rm = TRUE)

  # 7. Ratio été/hiver — indicateur persistant vs caduc
  ratio_summer_winter <- if (abs(ndvi_winter) > 0.01) ndvi_summer / ndvi_winter else NA_real_

  # 8. Variabilité intra-annuelle
  ndvi_cv <- sd(ndvi_ts, na.rm = TRUE) / mean(ndvi_ts, na.rm = TRUE)

  # 9. Asymétrie du cycle phénologique
  # Durée du verdissement vs durée de la sénescence
  if (!is.na(season$SOS) && !is.na(season$EOS) && !is.na(season$max_date)) {
    green_up_duration   <- season$max_date - season$SOS
    senescence_duration <- season$EOS - season$max_date
    asymmetry <- if (senescence_duration > 0) green_up_duration / senescence_duration else NA_real_
  } else {
    green_up_duration   <- NA_real_
    senescence_duration <- NA_real_
    asymmetry           <- NA_real_
  }

  # 10. Dérivées temporelles — nombre de pics
  # Détection de pics dans le profil NDVI
  diffs <- diff(sign(diff(ndvi_ts)))
  n_peaks <- sum(diffs == -2, na.rm = TRUE)

  # 11. Valeurs à des dates clés
  # DOY 90 (début avril), 135 (mi-mai), 182 (début juillet),
  # 244 (début sept), 305 (début nov)
  key_doys <- c(90, 135, 182, 244, 305)
  ndvi_at_keys <- sapply(key_doys, function(d) {
    idx <- which.min(abs(doy - d))
    ndvi_ts[idx]
  })

  # Assemblage du vecteur de features phénologiques
  metrics <- c(
    pheno_SOS               = season$SOS,
    pheno_EOS               = season$EOS,
    pheno_LOS               = season$LOS,
    pheno_max_doy           = season$max_date,
    pheno_max_ndvi          = ndvi_max,
    pheno_min_ndvi          = ndvi_min,
    pheno_amplitude         = amplitude,
    pheno_total_integral    = total_integral,
    pheno_growing_integral  = growing_integral,
    pheno_greening_rate     = greening_rate,
    pheno_browning_rate     = browning_rate,
    pheno_ndvi_winter       = ndvi_winter,
    pheno_ndvi_spring       = ndvi_spring,
    pheno_ndvi_summer       = ndvi_summer,
    pheno_ndvi_autumn       = ndvi_autumn,
    pheno_ratio_sum_win     = ratio_summer_winter,
    pheno_cv                = ndvi_cv,
    pheno_greenup_dur       = green_up_duration,
    pheno_senesc_dur        = senescence_duration,
    pheno_asymmetry         = asymmetry,
    pheno_n_peaks           = n_peaks,
    pheno_ndvi_doy090       = ndvi_at_keys[1],
    pheno_ndvi_doy135       = ndvi_at_keys[2],
    pheno_ndvi_doy182       = ndvi_at_keys[3],
    pheno_ndvi_doy244       = ndvi_at_keys[4],
    pheno_ndvi_doy305       = ndvi_at_keys[5]
  )

  metrics
}

# --- Extraction multi-indices ------------------------------------------------

#' Extraction de métriques phénologiques pour plusieurs indices spectraux
#' @param ts_data data.frame avec colonnes : date, et une colonne par indice
#' @param indices Noms des indices à traiter (défaut : NDVI, EVI)
#' @return Vecteur nommé de toutes les métriques
extract_multi_index_phenometrics <- function(ts_data, indices = c("NDVI", "EVI")) {
  all_metrics <- c()

  for (idx_name in indices) {
    if (!idx_name %in% names(ts_data)) next

    ts_values <- ts_data[[idx_name]]
    dates     <- ts_data$date

    # Lissage Savitzky-Golay
    ts_smooth <- smooth_savgol(ts_values,
                               order  = TS_PARAMS$sg_filter_order,
                               length = TS_PARAMS$sg_filter_length)

    # Extraction des métriques
    metrics <- extract_phenometrics(ts_smooth, dates)

    # Renommer avec le préfixe de l'indice
    names(metrics) <- gsub("^pheno_", paste0(idx_name, "_pheno_"), names(metrics))

    all_metrics <- c(all_metrics, metrics)
  }

  all_metrics
}

# --- Classification du type phénologique --------------------------------------

#' Détermination automatique du type phénologique (caduc / persistant)
#' @param ndvi_ts Série temporelle NDVI annuelle lissée
#' @param dates Dates correspondantes
#' @return Caractère : "deciduous", "evergreen", ou "semi_deciduous"
classify_phenotype <- function(ndvi_ts, dates) {
  if (all(is.na(ndvi_ts))) return(NA_character_)

  doy <- as.numeric(format(dates, "%j"))

  # NDVI hiver (DJF) vs été (JJA)
  winter <- doy <= 59 | doy >= 335
  summer <- doy >= 152 & doy <= 243

  ndvi_w <- mean(ndvi_ts[winter], na.rm = TRUE)
  ndvi_s <- mean(ndvi_ts[summer], na.rm = TRUE)

  # Coefficient de variation annuel
  cv <- sd(ndvi_ts, na.rm = TRUE) / mean(ndvi_ts, na.rm = TRUE)

  # Ratio été/hiver
  ratio <- if (ndvi_w > 0.05) ndvi_s / ndvi_w else Inf

  # Critères de classification
  if (cv < 0.10 && ratio < 1.3) {
    "evergreen"
  } else if (cv > 0.25 || ratio > 2.0) {
    "deciduous"
  } else {
    "semi_deciduous"
  }
}

# --- Décomposition en harmoniques de Fourier ----------------------------------

#' Décomposition harmonique d'une série temporelle NDVI
#' @param ndvi_ts Série temporelle NDVI
#' @param n_harmonics Nombre d'harmoniques à extraire (défaut 3)
#' @return Vecteur nommé avec amplitudes et phases des harmoniques
fourier_features <- function(ndvi_ts, n_harmonics = 3) {
  n <- length(ndvi_ts)
  if (n < 2 * n_harmonics || all(is.na(ndvi_ts))) {
    return(setNames(rep(NA_real_, 2 * n_harmonics + 1),
      c("fourier_mean", paste0("fourier_amp_", 1:n_harmonics),
        paste0("fourier_phase_", 1:n_harmonics))))
  }

  # FFT
  ft <- fft(ndvi_ts)

  # Amplitude et phase des premières harmoniques
  amplitudes <- Mod(ft[2:(n_harmonics + 1)]) * 2 / n
  phases     <- Arg(ft[2:(n_harmonics + 1)])

  features <- c(
    fourier_mean = mean(ndvi_ts, na.rm = TRUE),
    setNames(amplitudes, paste0("fourier_amp_", 1:n_harmonics)),
    setNames(phases, paste0("fourier_phase_", 1:n_harmonics))
  )

  features
}

