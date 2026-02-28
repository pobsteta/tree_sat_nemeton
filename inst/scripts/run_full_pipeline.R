#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI — Script complet : téléchargement, entraînement, prédiction
#
# Ce script :
#   1. Installe/charge le package treesatnemeton
#   2. Télécharge le dataset TreeSatAI-Time-Series depuis HuggingFace
#   3. Entraîne un modèle Random Forest (14 genres forestiers)
#   4. Produit une carte des essences sur votre AOI
#
# Modes d'utilisation :
#   - Données réelles TreeSatAI (défaut) : télécharge depuis HuggingFace
#   - Données synthétiques : USE_SYNTHETIC <- TRUE (test rapide, pas de téléchargement)
#
# Prérequis :
#   - R >= 4.1 avec les packages remotes, sf, terra, ranger installés
#   - Fichier data/aoi.gpkg pour la prédiction spatiale (étape 3)
#   - ~30 Go d'espace disque pour le dataset complet (mode réel)
#   - Connexion internet
# ==============================================================================

# --- Configuration -----------------------------------------------------------
USE_SYNTHETIC <- FALSE   # TRUE = test rapide avec données synthétiques
N_SAMPLES     <- 50      # Nombre d'échantillons par espèce (mode synthétique)
YEAR          <- 2021    # Année de simulation / prédiction
SKIP_PREDICT  <- FALSE   # TRUE = ne pas lancer la prédiction spatiale

# --- 0. Nettoyage et installation du package ---------------------------------
rm(list = ls(pattern = "^(?!USE_|N_SAMPLES|YEAR|SKIP_PREDICT)", perl = TRUE))
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

# Chemins locaux
PROJECT_ROOT <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "output")
MODELS_DIR   <- file.path(OUTPUT_DIR, "models")

# Créer l'arborescence du projet
init_project_dirs(PROJECT_ROOT)

# --- 1. Téléchargement du dataset (mode réel uniquement) ---------------------
if (!USE_SYNTHETIC) {
  cat("\n")
  cat("================================================================\n")
  cat("  Étape 1 — Téléchargement TreeSatAI depuis HuggingFace\n")
  cat("================================================================\n\n")

  # Labels + split + geojson (~10 Mo) — indispensable
  # sentinel-ts (~30 Go) — séries temporelles S1+S2 (optionnel, fallback synthétique)
  download_treesatai_hf(
    dest_dir   = file.path(DATA_DIR, "treesatai"),
    components = c("labels", "split", "geojson")
  )

  # Ajouter "sentinel-ts" aux components pour télécharger les séries temporelles :
  # download_treesatai_hf(
  #   dest_dir   = file.path(DATA_DIR, "treesatai"),
  #   components = c("labels", "split", "geojson", "sentinel-ts")
  # )

  cat("\nDonnées TreeSatAI prêtes.\n\n")
} else {
  cat("\n")
  cat("================================================================\n")
  cat("  Mode synthétique — pas de téléchargement\n")
  cat("================================================================\n\n")
}

# --- 2. Entraînement du modèle -----------------------------------------------
cat("================================================================\n")
cat("  Étape 2 — Entraînement du modèle Random Forest\n")
cat("================================================================\n\n")

# Supprimer l'ancien modèle s'il existe
old_model <- file.path(MODELS_DIR, "treesatai_rf.rds")
if (file.exists(old_model)) {
  file.remove(old_model)
  cat("Ancien modèle supprimé.\n\n")
}

# Lancer l'entraînement
# train_treesatai() gère automatiquement :
#   - Le chargement des labels et du split TreeSatAI
#   - La lecture des patches HDF5 ou la génération synthétique (fallback)
#   - Le filtrage des classes non-forestières (Cleared, etc.)
#   - La construction de la matrice de features
#   - L'entraînement, l'évaluation et les visualisations
if (USE_SYNTHETIC) {
  result_train <- train_treesatai(
    data_path  = NULL,
    mode       = "rf",
    n_samples  = N_SAMPLES,
    year       = YEAR,
    skip_viz   = FALSE
  )
} else {
  result_train <- train_treesatai(
    data_path  = file.path(DATA_DIR, "treesatai"),
    mode       = "rf",
    skip_viz   = FALSE
  )
}

# Résumé de l'entraînement
cat("\n")
cat("================================================================\n")
cat("  Résultats de l'entraînement\n")
cat("================================================================\n\n")

if (!is.null(result_train$evaluation)) {
  cat("  OA    :", round(result_train$evaluation$overall_accuracy * 100, 1), "%\n")
  cat("  Kappa :", round(result_train$evaluation$kappa, 3), "\n")
  if (!is.null(result_train$evaluation$macro_f1)) {
    cat("  F1    :", round(result_train$evaluation$macro_f1 * 100, 1), "%\n")
  }
  cat("  Features :", length(result_train$feature_cols), "\n")
  cat("  Modèle   :", file.path(MODELS_DIR, "treesatai_rf.rds"), "\n\n")
} else {
  cat("  Pas d'évaluation disponible.\n\n")
}

# --- 3. Prédiction spatiale sur une AOI (optionnel) --------------------------
if (!SKIP_PREDICT) {
  aoi_file <- file.path(DATA_DIR, "aoi.gpkg")

  if (file.exists(aoi_file)) {
    cat("================================================================\n")
    cat("  Étape 3 — Prédiction spatiale sur l'AOI\n")
    cat("================================================================\n\n")

    result <- predict_species_map(
      aoi_path      = aoi_file,
      auto_download = TRUE,
      year          = YEAR,
      resolution    = 10
    )

    cat("\n")
    cat("================================================================\n")
    cat("  Sorties de la prédiction\n")
    cat("================================================================\n\n")
    cat("  Carte des essences : output/carte_essences.tif\n")
    cat("  Carte confiance    : output/carte_confiance.tif\n")
    cat("  Polygones          : output/carte_essences.gpkg\n")
    cat("  Statistiques       : output/statistiques_essences.csv\n\n")
  } else {
    cat("================================================================\n")
    cat("  Étape 3 — Prédiction spatiale IGNORÉE\n")
    cat("================================================================\n\n")
    cat("  Fichier AOI non trouvé :", aoi_file, "\n")
    cat("  Pour lancer la prédiction, créez un fichier data/aoi.gpkg\n")
    cat("  avec le polygone de votre zone d'intérêt, puis relancez :\n\n")
    cat("    predict_species_map(\"data/aoi.gpkg\", auto_download = TRUE)\n\n")
  }
}

cat("================================================================\n")
cat("  Pipeline terminé !\n")
cat("================================================================\n")
