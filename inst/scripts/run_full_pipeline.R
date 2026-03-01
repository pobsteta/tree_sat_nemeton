#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI — Script complet : téléchargement, entraînement, prédiction
#
# Ce script :
#   1. Installe/charge le package treesatnemeton
#   2. Télécharge le dataset TreeSatAI-Time-Series depuis HuggingFace
#   3. Entraîne un modèle Random Forest (10 classes forestières)
#      — regroupement automatique des 20 espèces en 10 classes
#      — pondération des classes déséquilibrées
#      — sélection de features par Boruta
#      — support multi-année (features moyennées sur N ans)
#   4. Produit une carte des essences sur votre AOI
#
# Modes d'utilisation :
#   - Données réelles TreeSatAI (défaut) : télécharge depuis HuggingFace
#   - Données synthétiques : USE_SYNTHETIC <- TRUE (test rapide)
#
# Prérequis :
#   - R >= 4.1 avec les packages remotes, sf, terra, ranger installés
#   - Boruta recommandé : install.packages("Boruta")
#   - Fichier data/aoi.gpkg pour la prédiction spatiale (étape 3)
#   - ~30 Go d'espace disque pour le dataset complet (mode réel)
#   - Connexion internet
# ==============================================================================

# --- Configuration -----------------------------------------------------------
USE_SYNTHETIC <- FALSE   # TRUE = test rapide avec données synthétiques
N_SAMPLES     <- 50      # Nombre d'échantillons par espèce (mode synthétique)
YEARS         <- c(2019, 2020, 2021)  # Vecteur d'années (multi-année)
                                       # Mettre NULL ou c(2021) pour mono-année
USE_BORUTA    <- TRUE    # Sélection de features par Boruta
SKIP_PREDICT  <- FALSE   # TRUE = ne pas lancer la prédiction spatiale

# --- 0. Nettoyage et installation du package ---------------------------------
all_objs <- ls()
keep <- c("USE_SYNTHETIC", "N_SAMPLES", "YEARS", "USE_BORUTA", "SKIP_PREDICT")
rm(list = setdiff(all_objs, keep))
rm(all_objs, keep)
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

  download_treesatai_hf(
    dest_dir   = file.path(DATA_DIR, "treesatai"),
    components = c("labels", "split", "geojson")
  )

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
cat("================================================================\n")

use_multiyear <- !is.null(YEARS) && length(YEARS) > 1
if (use_multiyear) {
  cat("  Multi-année :", paste(YEARS, collapse = ", "), "\n")
}
if (USE_BORUTA) {
  cat("  Sélection de features : Boruta\n")
}
cat("  Classes : 10 groupes forestiers\n")
cat("  Pondération : inverse frequency weighting\n\n")

# Supprimer l'ancien modèle s'il existe
old_model <- file.path(MODELS_DIR, "treesatai_rf.rds")
if (file.exists(old_model)) {
  file.remove(old_model)
  cat("Ancien modèle supprimé.\n\n")
}

# Lancer l'entraînement
# train_treesatai() gère automatiquement :
#   - Le regroupement 20 espèces → 10 classes
#   - La pondération des classes déséquilibrées
#   - La sélection de features par Boruta (si use_boruta = TRUE)
#   - La construction multi-année (si years est un vecteur > 1 an)
#   - L'entraînement, l'évaluation et les visualisations
if (USE_SYNTHETIC) {
  result_train <- train_treesatai(
    data_path    = NULL,
    mode         = "rf",
    n_samples    = N_SAMPLES,
    years        = YEARS,
    use_boruta   = USE_BORUTA,
    skip_viz     = FALSE
  )
} else {
  result_train <- train_treesatai(
    data_path    = file.path(DATA_DIR, "treesatai"),
    mode         = "rf",
    use_boruta   = USE_BORUTA,
    skip_viz     = FALSE
  )
}

# Résumé de l'entraînement
cat("\n")
cat("================================================================\n")
cat("  Résultats de l'entraînement\n")
cat("================================================================\n\n")

if (!is.null(result_train$evaluation)) {
  cat("  OA       :", round(result_train$evaluation$overall_accuracy * 100, 1), "%\n")
  cat("  Kappa    :", round(result_train$evaluation$kappa, 3), "\n")
  if (!is.null(result_train$evaluation$macro_f1)) {
    cat("  F1       :", round(result_train$evaluation$macro_f1 * 100, 1), "%\n")
  }
  cat("  Features :", length(result_train$feature_cols), "\n")
  if (use_multiyear) {
    cat("  Années   :", paste(YEARS, collapse = ", "), "\n")
  }
  if (!is.null(result_train$boruta)) {
    cat("  Boruta   : actif\n")
  }
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
      year          = if (use_multiyear) max(YEARS) else YEARS[1],
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
