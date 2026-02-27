#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Script de test interactif pour RStudio
# Utilisation avec un fichier aoi.gpkg (zone d'intérêt)
#
# COMMENT UTILISER CE SCRIPT :
#   1. Ouvrir le projet dans RStudio (File > Open Project)
#   2. Modifier le chemin AOI_PATH ci-dessous
#   3. Exécuter bloc par bloc avec Ctrl+Enter
# ==============================================================================

# ==============================================================================
# BLOC 0 — CONFIGURATION (modifier ici)
# ==============================================================================

# Chemin vers votre fichier GeoPackage (zone d'intérêt)
AOI_PATH <- "~/mon_projet/aoi.gpkg"
# Ou en chemin absolu :
# AOI_PATH <- "/home/utilisateur/donnees/aoi.gpkg"

# Année d'analyse
YEAR <- 2021

# Nombre d'échantillons par espèce (mode synthétique)
N_SAMPLES <- 30  # Réduire pour un test rapide, 50+ pour un vrai run

# Méthode : "rf" (Random Forest) ou "both" (RF + CNN)
METHOD <- "rf"

# ==============================================================================
# BLOC 1 — Installation des packages
# ==============================================================================

# Exécuter une seule fois — installe les packages manquants
install_if_needed <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installation de : ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cran.r-project.org")
  } else {
    message("Tous les packages sont déjà installés.")
  }
}

# Packages essentiels (sans torch pour commencer)
install_if_needed(c(
  "sf", "terra", "here", "cli", "glue",
  "dplyr", "tidyr", "purrr", "readr", "stringr", "lubridate",
  "ggplot2", "patchwork", "viridis", "RColorBrewer", "pheatmap",
  "ranger", "signal", "pracma", "jsonlite",
  "scales", "data.table"
))

# ==============================================================================
# BLOC 2 — Chargement des modules du projet
# ==============================================================================

# Définir la racine du projet (là où est le dossier R/)
# Option A : si vous avez ouvert le .Rproj ou êtes dans le bon dossier
setwd(here::here())

# Option B : chemin absolu si nécessaire
# setwd("/chemin/vers/tree_sat_nemeton")

# Charger tous les modules
source("R/00_config.R")
source("R/01_utils.R")
source("R/02_phenology.R")
source("R/03_data_acquisition.R")
source("R/04_visualization.R")
source("R/05_classification.R")

message("\n=== Modules TreeSatAI-TS chargés ===\n")

# ==============================================================================
# BLOC 3 — Charger et explorer l'AOI
# ==============================================================================

# Charger le GeoPackage
aoi <- sf::st_read(AOI_PATH)

# Explorer la structure
cat("--- Structure de l'AOI ---\n")
print(str(aoi))

cat("\n--- Aperçu ---\n")
print(head(aoi))

cat("\n--- CRS ---\n")
print(sf::st_crs(aoi))

cat("\n--- Bbox ---\n")
print(sf::st_bbox(aoi))

cat("\n--- Géométries ---\n")
print(table(sf::st_geometry_type(aoi)))

cat("\n--- Colonnes disponibles ---\n")
print(names(aoi))

cat("\n--- Nombre d'entités ---\n")
print(nrow(aoi))

# Visualisation rapide de l'AOI
plot(sf::st_geometry(aoi), main = "Zone d'intérêt (AOI)", col = "lightgreen", border = "darkgreen")

# ==============================================================================
# BLOC 4 — Adapter l'AOI pour le pipeline
# ==============================================================================

# Le pipeline attend un champ "essence" (ou "species") et un champ "idp" (id parcelle)
# Adaptez les noms de colonnes ci-dessous à votre GeoPackage

cat("Colonnes de votre AOI :\n")
print(names(aoi))

# --- OPTION A : Votre AOI contient déjà des espèces identifiées ---
# Décommentez et adaptez selon vos noms de colonnes :

# aoi_prepared <- aoi |>
#   dplyr::rename(
#     idp          = NOM_DE_VOTRE_COLONNE_ID,      # ex: "id", "fid", "gid"
#     species_name = NOM_DE_VOTRE_COLONNE_ESPECE    # ex: "essence", "species", "espar"
#   ) |>
#   dplyr::mutate(idp = as.character(idp))

# --- OPTION B : Votre AOI est une zone sans espèces (juste une emprise) ---
# On génère des données synthétiques à l'intérieur de l'AOI

# --- OPTION C : Votre AOI contient des codes IFN ---
# Mapping des codes ESPAR de l'IFN vers les noms TreeSatAI
IFN_TO_TREESATAI <- c(
  "01" = "Chêne pédonculé",     "02" = "Chêne sessile",
  "03" = "Chêne pubescent",     "04" = "Chêne vert",
  "09" = "Hêtre",               "12" = "Châtaignier",
  "13" = "Charme",              "10" = "Bouleau verruqueux",
  "14" = "Frêne commun",        "16" = "Érable sycomore",
  "17" = "Peupliers",           "20" = "Robinier faux-acacia",
  "52" = "Épicéa commun",       "61" = "Sapin pectiné",
  "53" = "Douglas",             "51" = "Pin sylvestre",
  "54" = "Pin maritime",        "55" = "Pin noir",
  "56" = "Pin d'Alep",          "63" = "Mélèze d'Europe"
)

# ==============================================================================
# BLOC 5A — Test rapide avec DONNÉES SYNTHÉTIQUES (recommandé pour commencer)
# ==============================================================================

# Ce mode ne nécessite aucune image Sentinel-2
# Il génère des profils phénologiques réalistes pour tester tout le pipeline

cat("\n=== MODE SYNTHÉTIQUE ===\n")
cat("Génération de données de test...\n")

# Mettre à jour l'année dans la config
TS_PARAMS$year       <- YEAR
TS_PARAMS$start_date <- paste0(YEAR, "-01-01")
TS_PARAMS$end_date   <- paste0(YEAR, "-12-31")

# Générer le dataset synthétique
ts_long <- generate_synthetic_dataset(
  n_samples_per_species = N_SAMPLES,
  year = YEAR
)

# Ajouter les indices spectraux
ts_long <- ts_long |>
  dplyr::mutate(
    NDVI   = calc_ndvi(B08, B04),
    EVI    = calc_evi(B08, B04, B02),
    NDWI   = calc_ndwi(B08, B11),
    CRI    = calc_cri(B03, B05),
    RENDVI = calc_rendvi(B08, B05),
    NBR    = calc_nbr(B08, B12)
  )

cat(sprintf("Dataset : %d lignes, %d parcelles, %d espèces\n",
            nrow(ts_long),
            length(unique(ts_long$plot_id)),
            length(unique(ts_long$species_name))))

# Aperçu rapide
print(head(ts_long[, c("plot_id", "species_name", "date", "NDVI", "EVI")]))

# ==============================================================================
# BLOC 5B — Mode DONNÉES RÉELLES (si vous avez des images Sentinel-2)
# ==============================================================================

# Décommentez ce bloc si vous avez des images S2 et des parcelles avec espèces

# # Charger les parcelles depuis l'AOI (avec champ espèce)
# plots <- aoi_prepared  # Depuis le BLOC 4
#
# # Répertoire contenant les images Sentinel-2 L2A
# S2_DIR <- "~/donnees/sentinel2/L2A/"
#
# # Extraction des séries temporelles
# ts_long <- extract_s2_timeseries(
#   plots   = plots,
#   s2_dir  = S2_DIR,
#   bands   = S2_BAND_NAMES,
#   buffer_m = 15   # rayon du buffer (m) autour de chaque point
# )
#
# # Ajouter les indices
# ts_long <- ts_long |>
#   dplyr::mutate(
#     NDVI   = calc_ndvi(B08, B04),
#     EVI    = calc_evi(B08, B04, B02),
#     NDWI   = calc_ndwi(B08, B11),
#     CRI    = calc_cri(B03, B05),
#     RENDVI = calc_rendvi(B08, B05),
#     NBR    = calc_nbr(B08, B12)
#   )

# ==============================================================================
# BLOC 6 — Visualisation des profils phénologiques
# ==============================================================================

cat("\n=== PROFILS PHÉNOLOGIQUES ===\n")

# Profils NDVI de toutes les espèces
p1 <- plot_phenological_profiles(ts_long, index_name = "NDVI")
print(p1)

# Comparaison caducifolié vs sempervirent
p2 <- plot_deciduous_vs_evergreen(ts_long)
print(p2)

# Profils par type (feuillu/résineux × caduc/persistant)
p3 <- plot_profiles_by_type(ts_long)
print(p3)

# --- Zoom sur quelques espèces intéressantes ---
species_focus <- c("Hêtre", "Épicéa commun", "Mélèze d'Europe", "Chêne vert")

ts_focus <- ts_long |>
  dplyr::filter(species_name %in% species_focus)

p_focus <- plot_phenological_profiles(ts_focus, index_name = "NDVI")
print(p_focus + ggplot2::ggtitle("4 profils contrastés : caduc, persistant, mélèze, chêne vert"))

# ==============================================================================
# BLOC 7 — Construction de la matrice de features
# ==============================================================================

cat("\n=== EXTRACTION DES FEATURES ===\n")

feature_matrix <- build_feature_matrix(ts_long)

cat(sprintf("Matrice de features : %d parcelles × %d colonnes\n",
            nrow(feature_matrix), ncol(feature_matrix)))

# Aperçu des types de features
feature_types <- data.frame(
  type = c("Bandes temporelles", "Indices temporels",
           "Stats bandes", "Stats indices", "Phénologie", "Fourier"),
  pattern = c("^B\\d{2}_d|^B8A_d", "^(NDVI|EVI|NDWI|CRI|RENDVI|NBR)_d",
              "^B\\d{2}_(mean|sd)", "^(NDVI|EVI)_(mean|sd)",
              "pheno_", "fourier_")
)

for (i in seq_len(nrow(feature_types))) {
  n <- sum(grepl(feature_types$pattern[i], names(feature_matrix)))
  cat(sprintf("  %-25s : %d features\n", feature_types$type[i], n))
}

# ==============================================================================
# BLOC 8 — Classification Random Forest
# ==============================================================================

cat("\n=== CLASSIFICATION RANDOM FOREST ===\n")

# Séparation train/test
split <- split_train_test(feature_matrix, target_col = "species_name")
train_data <- split$train
test_data  <- split$test

# Sélection des features
feature_cols <- select_features(feature_matrix)
cat(sprintf("Features sélectionnées : %d\n", length(feature_cols)))

# Entraînement
rf_model <- train_random_forest(train_data, feature_cols)

# Prédiction
rf_preds <- predict_rf(rf_model, test_data)

# Évaluation
rf_eval <- evaluate_classification(
  y_true      = test_data$species_name,
  y_pred      = rf_preds$predicted_class,
  class_names = SPECIES$french
)

# ==============================================================================
# BLOC 9 — Visualisations des résultats
# ==============================================================================

cat("\n=== RÉSULTATS ===\n")

# Matrice de confusion
p_conf <- plot_confusion_matrix(
  rf_eval$confusion_matrix,
  class_names = SPECIES$french,
  title = "Matrice de confusion — 20 espèces"
)
print(p_conf)

# Importance des variables (top 30)
rf_importance <- get_variable_importance(rf_model)
p_imp <- plot_variable_importance(rf_importance, top_n = 30)
print(p_imp)

# Métriques par espèce
p_met <- plot_species_metrics(rf_eval$per_class)
print(p_met)

# Heatmap des signatures NDVI
plot_spectral_heatmap(feature_matrix)

# ==============================================================================
# BLOC 10 — Sauvegarder les résultats
# ==============================================================================

cat("\n=== SAUVEGARDE ===\n")

# Créer les dossiers de sortie
dir.create("output", showWarnings = FALSE)
dir.create("output/models", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Modèle
save_model(rf_model, rf_eval, model_name = "treesatai_rf_test")

# Figures
ggsave("figures/01_phenological_profiles.png", p1, width = 20, height = 15, units = "cm", dpi = 300)
ggsave("figures/02_deciduous_vs_evergreen.png", p2, width = 20, height = 12, units = "cm", dpi = 300)
ggsave("figures/04_confusion_matrix.png", p_conf, width = 25, height = 22, units = "cm", dpi = 300)
ggsave("figures/05_variable_importance.png", p_imp, width = 20, height = 20, units = "cm", dpi = 300)
ggsave("figures/06_species_metrics.png", p_met, width = 20, height = 12, units = "cm", dpi = 300)

# CSV
readr::write_csv(rf_eval$per_class, "output/per_class_metrics.csv")
readr::write_csv(rf_importance, "output/variable_importance.csv")
readr::write_csv(feature_matrix, "output/feature_matrix.csv")

cat("\nFichiers sauvegardés dans output/ et figures/\n")

# ==============================================================================
# BLOC 11 — Exploration interactive
# ==============================================================================

cat("\n=== EXPLORATION INTERACTIVE ===\n")
cat("Objets disponibles dans l'environnement :\n")
cat("  aoi            — votre zone d'intérêt (sf)\n")
cat("  ts_long        — séries temporelles longues\n")
cat("  feature_matrix — matrice de features (1 ligne/parcelle)\n")
cat("  rf_model       — modèle Random Forest entraîné\n")
cat("  rf_eval        — résultats d'évaluation\n")
cat("  rf_importance  — importance des variables\n")
cat("  SPECIES        — tableau des 20 espèces\n")
cat("\nExemples de commandes à tester :\n")
cat('  View(SPECIES)                           # Voir les 20 espèces\n')
cat('  View(rf_eval$per_class)                 # Métriques par espèce\n')
cat('  head(rf_importance, 20)                 # Top 20 features\n')
cat('  table(feature_matrix$species_name)      # Distribution des classes\n')

# ==============================================================================
# BLOC 12 (optionnel) — Prédire sur de nouvelles parcelles
# ==============================================================================

# Si vous avez de nouvelles parcelles sans étiquette d'espèce,
# vous pouvez utiliser le modèle entraîné pour prédire :

# new_parcels <- sf::st_read("chemin/vers/nouvelles_parcelles.gpkg")
# new_ts <- extract_s2_timeseries(new_parcels, S2_DIR)
# new_features <- build_feature_matrix(new_ts)
# predictions <- predict_rf(rf_model, new_features)
#
# # Résultat : espèce prédite + probabilités
# new_parcels$espece_predite <- predictions$predicted_class
# new_parcels$proba_max <- apply(predictions$probabilities, 1, max)
#
# # Visualiser
# plot(new_parcels["espece_predite"], main = "Espèces prédites")
#
# # Exporter
# sf::st_write(new_parcels, "output/predictions.gpkg")
