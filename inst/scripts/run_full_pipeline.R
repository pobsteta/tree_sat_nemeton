# ==============================================================================
# TreeSatAI — Script complet : téléchargement, entraînement, prédiction
#
# Ce script :
#   1. Installe/charge le package treesatnemeton
#   2. Télécharge le dataset TreeSatAI-Time-Series depuis HuggingFace
#   3. Entraîne un modèle Random Forest sur les données réelles
#   4. Produit une carte des essences sur votre AOI
#
# Prérequis :
#   - R >= 4.1 avec les packages remotes, sf, terra, ranger installés
#   - Fichier data/aoi.gpkg dans votre répertoire de travail
#   - ~30 Go d'espace disque pour le dataset complet
#   - Connexion internet
# ==============================================================================

# --- 0. Nettoyage et installation du package ---------------------------------
rm(list = ls())
gc()

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Installer la dernière version du package
remotes::install_github(
  "pobsteta/tree_sat_nemeton",
  ref = "claude/treesatai-species-classification-HB6Vr",
  force = TRUE,
  upgrade = "never"
)

library(treesatnemeton)

# --- 1. Télécharger le dataset TreeSatAI depuis HuggingFace -----------------
# Labels + split + geojson (~10 Mo) — indispensable
# sentinel-ts (~30 Go) — séries temporelles S1+S2 sur 1 an complet
download_treesatai_hf(
  components = c("labels", "split", "geojson", "sentinel-ts")
)

# Vérifier les labels
labels <- load_treesatai_labels()
split  <- load_treesatai_split()

cat("\n=== Dataset TreeSatAI ===\n")
cat("Patchs labelisés :", length(unique(labels$patch_id)), "\n")
cat("Genres uniques   :", length(unique(labels$genus)), "\n")
cat("Train            :", length(split$train), "patchs\n")
cat("Test             :", length(split$test), "patchs\n\n")

# --- 2. Entraîner un modèle sur les données réelles -------------------------
# Supprime l'ancien modèle synthétique s'il existe
old_model <- file.path(MODELS_DIR, "treesatai_rf.rds")
if (file.exists(old_model)) {
  file.remove(old_model)
  cat("Ancien modèle synthétique supprimé.\n")
}

# Entraînement sur les données TreeSatAI réelles
result_train <- train_treesatai(
  data_path = file.path(DATA_DIR, "treesatai"),
  mode      = "rf",
  skip_viz  = FALSE  # Générer les visualisations (profils, confusion, etc.)
)

cat("\n=== Modèle entraîné ===\n")
cat("OA    :", round(result_train$evaluation$overall_accuracy * 100, 1), "%\n")
cat("Kappa :", round(result_train$evaluation$kappa, 3), "\n")
cat("F1    :", round(result_train$evaluation$macro_f1 * 100, 1), "%\n\n")

# --- 3. Prédire sur votre AOI -----------------------------------------------
result <- predict_species_map(
  aoi_path      = "data/aoi.gpkg",
  auto_download = TRUE,
  year          = 2023,
  resolution    = 10
)

cat("\n=== Terminé ===\n")
cat("Carte des essences : output/carte_essences.tif\n")
cat("Carte confiance    : output/carte_confiance.tif\n")
cat("Polygones          : output/carte_essences.gpkg\n")
cat("Statistiques       : output/statistiques_essences.csv\n")
