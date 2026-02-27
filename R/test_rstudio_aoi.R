#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Script interactif pour RStudio
#
# OBJECTIF : Donner un fichier aoi.gpkg → obtenir la carte des essences
#
# COMMENT UTILISER :
#   1. Ouvrir le projet dans RStudio (File > Open Project)
#   2. Modifier AOI_PATH ci-dessous
#   3. Configurer vos identifiants CDSE (voir ci-dessous)
#   4. Exécuter bloc par bloc avec Ctrl+Enter (ou Ctrl+Shift+Enter par bloc)
#
# MODES DISPONIBLES :
#   - "download" (défaut) → télécharge S2 (+S1) depuis Copernicus CDSE
#   - "local"             → images Sentinel déjà téléchargées sur disque
#   - "pytorch"           → utilise un modèle Deep Learning (TempCNN, etc.)
# ==============================================================================

# ==============================================================================
# BLOC 0 — CONFIGURATION (seule partie à modifier)
# ==============================================================================

# Chemin vers votre fichier GeoPackage (zone d'intérêt)
AOI_PATH <- "~/mon_projet/aoi.gpkg"

# ----- MODE DE FONCTIONNEMENT (décommenter UN seul mode) -----

# MODE 1 (défaut) : Téléchargement automatique S2 (+S1) depuis Copernicus
#   Nécessite un compte gratuit sur https://dataspace.copernicus.eu
#   Configurer dans ~/.Renviron :
#     CDSE_USERNAME=votre@email.com
#     CDSE_PASSWORD=motdepasse
MODE <- "download"

# MODE 2 : Données locales (vous avez déjà les images)
# MODE <- "local"

# MODE 3 : PyTorch (utilise un modèle DL au lieu de Random Forest)
#   Nécessite : conda activate treesat (ou setup_python_env() depuis R)
# MODE <- "pytorch"

# ----- PARAMÈTRES -----

# (Mode local) Répertoire de vos images Sentinel-2 L2A
S2_DIR <- NULL
# S2_DIR <- "~/donnees/sentinel2/L2A/"

# (Mode local) Répertoire de vos images Sentinel-1 GRD
S1_DIR <- NULL
# S1_DIR <- "~/donnees/sentinel1/GRD/"

# Inclure Sentinel-1 (radar) ? Utile si beaucoup de nuages
USE_S1 <- FALSE

# (Optionnel) Modèle pré-entraîné
MODEL_PATH <- NULL
# MODEL_PATH <- "output/models/treesatai_rf.rds"

# Année d'analyse
YEAR <- 2023

# Résolution cible (mètres) — 10m ou 20m
RESOLUTION <- 10

# ==============================================================================
# BLOC 1 — Installation des packages (une seule fois)
# ==============================================================================

install_if_needed <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installation de : ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cran.r-project.org")
  } else {
    message("Tous les packages sont installés.")
  }
}

install_if_needed(c(
  "sf", "terra", "here", "cli", "glue", "httr2", "rstac",
  "dplyr", "tidyr", "purrr", "readr", "stringr", "lubridate",
  "ggplot2", "patchwork", "viridis", "RColorBrewer", "pheatmap",
  "ranger", "signal", "pracma", "jsonlite", "data.table", "scales"
))

# ==============================================================================
# BLOC 2 — Chargement du projet
# ==============================================================================

setwd(here::here())  # Se placer à la racine du projet

# Charger tous les modules
source("R/00_config.R")
source("R/01_utils.R")
source("R/02_phenology.R")
source("R/03_data_acquisition.R")
source("R/04_visualization.R")
source("R/05_classification.R")
source("R/07_predict_aoi.R")

cat("\n=== TreeSatAI-TS chargé — prêt à cartographier ===\n\n")

# ==============================================================================
# BLOC 3 — Explorer l'AOI
# ==============================================================================

aoi <- sf::st_read(AOI_PATH)

cat("--- Zone d'intérêt ---\n")
cat("Fichier     :", AOI_PATH, "\n")
cat("Entités     :", nrow(aoi), "\n")
cat("Géométrie   :", as.character(unique(sf::st_geometry_type(aoi))), "\n")
cat("CRS         :", sf::st_crs(aoi)$input, "\n")
cat("Colonnes    :", paste(names(aoi), collapse = ", "), "\n")

# Surface
aoi_area <- sf::st_area(sf::st_transform(aoi, 2154))
cat("Surface     :", round(sum(as.numeric(aoi_area)) / 10000, 1), "ha\n")

# Estimation du temps
n_pixels_est <- round(sum(as.numeric(aoi_area)) / RESOLUTION^2)
cat("Pixels est. :", format(n_pixels_est, big.mark = " "), paste0("(", RESOLUTION, "m)\n"))

# Visualisation
plot(sf::st_geometry(aoi), main = "Zone d'intérêt", col = "lightgreen", border = "darkgreen")

# ==============================================================================
# BLOC 4 — LANCER LA DÉTECTION DES ESSENCES
# ==============================================================================

# La commande s'adapte automatiquement au MODE choisi

if (MODE == "download") {
  # --- Télécharger S2 (+S1) puis classifier ---
  result <- predict_species_map(
    aoi_path      = AOI_PATH,
    model_path    = MODEL_PATH,
    year          = YEAR,
    output_dir    = file.path(here::here(), "output"),
    resolution    = RESOLUTION,
    auto_download = TRUE,
    use_s1        = USE_S1
  )

} else if (MODE == "local") {
  # --- Données locales ---
  result <- predict_species_map(
    aoi_path   = AOI_PATH,
    s2_dir     = S2_DIR,
    s1_dir     = S1_DIR,
    model_path = MODEL_PATH,
    year       = YEAR,
    output_dir = file.path(here::here(), "output"),
    resolution = RESOLUTION,
    use_s1     = USE_S1
  )

} else if (MODE == "pytorch") {
  # --- Mode PyTorch (Deep Learning) ---
  result <- predict_species_map(
    aoi_path      = AOI_PATH,
    s2_dir        = S2_DIR,
    model_path    = MODEL_PATH,
    year          = YEAR,
    output_dir    = file.path(here::here(), "output"),
    resolution    = RESOLUTION,
    auto_download = FALSE,
    use_pytorch   = TRUE,
    pytorch_model = "tempcnn"
  )

} else {
  stop(paste0(
    "Mode '", MODE, "' non reconnu.\n",
    "Modes valides : 'download', 'local', 'pytorch'"
  ))
}

# ==============================================================================
# BLOC 5 — Visualiser les résultats
# ==============================================================================

# --- 5a. Carte des essences ---
terra::plot(result$species_raster,
            main = "Carte des essences forestières",
            col = SPECIES_COLORS)

# --- 5b. Carte de confiance ---
terra::plot(result$confidence_raster,
            main = "Confiance de la prédiction",
            col = viridis::viridis(100))

# --- 5c. Tableau récapitulatif ---
cat("\n=== Essences détectées ===\n\n")
detected <- result$statistics[result$statistics$n_pixels > 0, ]
detected_display <- detected[, c("espece", "surface_ha", "pct", "confiance_moy", "type")]
names(detected_display) <- c("Espèce", "Surface (ha)", "%", "Confiance (%)", "Type")
print(detected_display, row.names = FALSE)

# --- 5d. Graphique en barres ---
library(ggplot2)

p_barres <- ggplot(detected, aes(x = reorder(espece, -surface_ha),
                                  y = surface_ha, fill = type)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = c("feuillu" = "#4daf4a", "résineux" = "#377eb8"),
                    name = "Type") +
  labs(title = "Surface par espèce détectée",
       x = NULL, y = "Surface (ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
print(p_barres)

# ==============================================================================
# BLOC 6 — Exporter les résultats
# ==============================================================================

cat("\n=== Fichiers de sortie ===\n")
cat("  output/carte_essences.tif    — raster classifié (ouvrir dans QGIS)\n")
cat("  output/carte_confiance.tif   — carte de confiance\n")
cat("  output/carte_essences.gpkg   — polygones par espèce (ouvrir dans QGIS)\n")
cat("  output/statistiques_essences.csv — statistiques\n")
cat("  output/legende_especes.csv   — légende des codes\n")

# Pour ouvrir dans QGIS directement depuis R :
# system("qgis output/carte_essences.gpkg &")

# ==============================================================================
# BLOC 7 (optionnel) — Télécharger les données satellite
# ==============================================================================

# Si vous voulez UNIQUEMENT télécharger les données (sans classifier) :

# source("R/08_download_satellite.R")
#
# # Télécharger Sentinel-2 + Sentinel-1 sur votre AOI
# sat_data <- download_satellite_data(
#   aoi_path    = AOI_PATH,
#   year        = YEAR,
#   download_s2 = TRUE,    # Sentinel-2 optique
#   download_s1 = TRUE     # Sentinel-1 radar
# )
#
# # Ou séparément :
# s2_dir <- download_s2_for_aoi(aoi, year = YEAR, max_cloud = 20)
# s1_dir <- download_s1_for_aoi(aoi, year = YEAR)
#
# # Puis classifier avec les données téléchargées :
# result <- predict_species_map(
#   aoi_path = AOI_PATH,
#   s2_dir   = file.path(sat_data$s2_dir, "bands"),
#   s1_dir   = sat_data$s1_dir,
#   use_s1   = TRUE,
#   year     = YEAR
# )

# ==============================================================================
# BLOC 7b (optionnel) — Utiliser un modèle PyTorch (TempCNN, Transformer...)
# ==============================================================================

# ÉTAPE 1 : Installer l'environnement Python conda "treesat"
#   source("R/09_python_bridge.R")
#   setup_python_env()   # crée l'env conda depuis environment.yml

# ÉTAPE 2 : Entraîner un modèle PyTorch
#   # Depuis R :
#   py_train_model("data/processed/feature_matrix.csv", model_type = "tempcnn")
#   # Ou depuis le terminal :
#   # conda activate treesat
#   # python python/train.py --data data/processed/feature_matrix.csv --model tempcnn

# ÉTAPE 3 : Prédire avec le modèle PyTorch
#   result <- predict_species_map(
#     aoi_path      = AOI_PATH,
#     auto_download = TRUE,
#     year          = 2023,
#     use_pytorch   = TRUE,
#     pytorch_model = "tempcnn",   # ou "lstm", "transformer", "inception"
#     model_path    = "output/models/treesatai_tempcnn_best.pt"
#   )

# Architectures disponibles :
#   - "tempcnn"     : Temporal CNN (Pelletier et al. 2019) — rapide, bon par défaut
#   - "lstm"        : LSTM bidirectionnel — capture les dépendances longues
#   - "transformer" : Transformer encoder — attention sur les dates clés
#   - "inception"   : InceptionTime — multi-échelle temporelle
#   - "multisource" : TempCNN fusion S2+S1 (nécessite use_s1 = TRUE)

# ==============================================================================
# BLOC 8 (optionnel) — Filtrage par confiance
# ==============================================================================

# Garder uniquement les pixels avec une confiance > 60%
if (!is.null(result)) {
  seuil_confiance <- 0.6

  carte_filtree <- result$species_raster
  carte_filtree[result$confidence_raster < seuil_confiance] <- NA

  terra::plot(carte_filtree,
              main = paste0("Essences (confiance > ", seuil_confiance * 100, "%)"),
              col = SPECIES_COLORS)

  # Pourcentage de pixels conservés
  n_total <- sum(!is.na(terra::values(result$species_raster)))
  n_filtre <- sum(!is.na(terra::values(carte_filtree)))
  cat(sprintf("\nPixels conservés : %d / %d (%.1f%%)\n",
              n_filtre, n_total, n_filtre / n_total * 100))
}

# ==============================================================================
# BLOC 9 (optionnel) — Superposer AOI + carte dans un même plot
# ==============================================================================

if (!is.null(result)) {
  terra::plot(result$species_raster,
              main = "Essences détectées dans la zone d'intérêt",
              col = SPECIES_COLORS)
  plot(sf::st_geometry(sf::st_transform(aoi, terra::crs(result$species_raster))),
       add = TRUE, border = "red", lwd = 2)
  legend("bottomright", legend = "AOI", col = "red", lty = 1, lwd = 2, cex = 0.8)
}
