#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series (IGNF) — Pipeline CLI
# Classification de 10 classes forestières européennes
# par séries temporelles Sentinel-1/2 (mono ou multi-année)
#
# Usage :
#   Rscript inst/scripts/06_pipeline.R --synthetic
#   Rscript inst/scripts/06_pipeline.R --synthetic --years 2019,2020,2021
#   Rscript inst/scripts/06_pipeline.R --synthetic --no-boruta
#   Rscript inst/scripts/06_pipeline.R --data /chemin/dataset
#   Rscript inst/scripts/06_pipeline.R --data /chemin --mode cnn
# ==============================================================================

# --- Chargement des modules ---------------------------------------------------
source(file.path(here::here(), "R", "00_config.R"))
source(file.path(here::here(), "R", "01_utils.R"))
source(file.path(here::here(), "R", "02_phenology.R"))
source(file.path(here::here(), "R", "03_data_acquisition.R"))
source(file.path(here::here(), "R", "04_visualization.R"))
source(file.path(here::here(), "R", "05_classification.R"))

# --- Arguments en ligne de commande ------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

DATA_PATH     <- NULL
MODE          <- "rf"       # "rf", "cnn" ou "both"
SKIP_VIZ      <- FALSE
USE_SYNTHETIC <- FALSE
N_SAMPLES     <- 50         # Échantillons par espèce (mode synthétique)
YEARS         <- NULL        # NULL = single year, sinon vecteur d'années
USE_BORUTA    <- TRUE
BORUTA_MAXRUNS <- 100

for (i in seq_along(args)) {
  if (args[i] == "--data" && i < length(args))         DATA_PATH <- args[i + 1]
  if (args[i] == "--mode" && i < length(args))         MODE      <- args[i + 1]
  if (args[i] == "--no-viz")                            SKIP_VIZ  <- TRUE
  if (args[i] == "--synthetic")                         USE_SYNTHETIC <- TRUE
  if (args[i] == "--n-samples" && i < length(args))    N_SAMPLES <- as.integer(args[i + 1])
  if (args[i] == "--years" && i < length(args))        YEARS <- as.integer(strsplit(args[i + 1], ",")[[1]])
  if (args[i] == "--no-boruta")                         USE_BORUTA <- FALSE
  if (args[i] == "--boruta-maxruns" && i < length(args)) BORUTA_MAXRUNS <- as.integer(args[i + 1])
}

# ==============================================================================
cli::cli_h1("Pipeline TreeSatAI-Time-Series")
cli::cli_text("")
# ==============================================================================

# Résoudre le mode multi-année
use_multiyear <- !is.null(YEARS) && length(YEARS) > 1
if (is.null(YEARS)) YEARS <- 2021

# Nombre de classes (10 groupes)
n_classes <- nrow(SPECIES_GROUPS_INFO)
eval_class_names <- SPECIES_GROUPS_INFO$group

# ==============================================================================
# ÉTAPE 1 : ACQUISITION DES DONNÉES
# ==============================================================================
cli::cli_h2("Étape 1 — Acquisition des données")

if (USE_SYNTHETIC) {
  # ---- Mode données synthétiques TreeSatAI ----
  cli::cli_alert_info("Mode synthétique : génération de {N_SAMPLES} échantillons × {nrow(SPECIES)} espèces")
  if (use_multiyear) {
    cli::cli_alert_info("Multi-année : {length(YEARS)} ans ({paste(YEARS, collapse=', ')})")
  }
  cli::cli_text("")
  cli::cli_text("Les profils phénologiques sont basés sur les signatures spectrales")
  cli::cli_text("caractéristiques de chaque essence (double logistique + bruit).")
  cli::cli_text("")

  # Générer la première année pour la visualisation
  ts_long <- generate_synthetic_dataset(n_samples_per_species = N_SAMPLES, year = YEARS[1])
  cli::cli_alert_success("Dataset synthétique généré : {nrow(ts_long)} observations")

} else if (!is.null(DATA_PATH)) {
  # ---- Mode données réelles ----
  cli::cli_alert_info("Chargement des données depuis : {DATA_PATH}")

  if (dir.exists(DATA_PATH)) {
    plot_files <- list.files(DATA_PATH, pattern = "\\.(gpkg|shp|geojson)$",
                             recursive = TRUE, full.names = TRUE)
    if (length(plot_files) > 0) {
      plots <- sf::st_read(plot_files[1], quiet = TRUE)
      cli::cli_alert_success("{nrow(plots)} parcelles chargées")
    }

    ts_files <- list.files(DATA_PATH, pattern = "\\.(csv|parquet)$",
                           recursive = TRUE, full.names = TRUE)
    if (length(ts_files) > 0) {
      ts_long <- readr::read_csv(ts_files[1], show_col_types = FALSE)
      cli::cli_alert_success("Séries temporelles chargées : {nrow(ts_long)} lignes")
    } else {
      s2_dir <- list.dirs(DATA_PATH, recursive = FALSE)
      s2_dir <- s2_dir[grep("sentinel|S2", s2_dir, ignore.case = TRUE)]
      if (length(s2_dir) > 0) {
        ts_long <- extract_s2_timeseries(plots, s2_dir[1])
      } else {
        cli::cli_alert_danger("Aucune donnée Sentinel-2 trouvée dans {DATA_PATH}")
        stop("Données manquantes")
      }
    }
  } else {
    cli::cli_alert_danger("Répertoire introuvable : {DATA_PATH}")
    stop("Chemin invalide")
  }

} else {
  cli::cli_alert_danger("Aucun répertoire de données spécifié.")
  cli::cli_text("")
  cli::cli_text("Usage :")
  cli::cli_text("  Rscript inst/scripts/06_pipeline.R --synthetic")
  cli::cli_text("  Rscript inst/scripts/06_pipeline.R --synthetic --years 2019,2020,2021")
  cli::cli_text("  Rscript inst/scripts/06_pipeline.R --synthetic --no-boruta")
  cli::cli_text("  Rscript inst/scripts/06_pipeline.R --data /chemin/vers/donnees")
  stop("Argument --data ou --synthetic obligatoire.")
}

# Regroupement des 20 espèces en 10 classes forestières
cli::cli_alert_info("Regroupement des espèces en {n_classes} classes")
ts_long <- remap_species_groups(ts_long)

# ==============================================================================
# ÉTAPE 2 : VISUALISATION DES PROFILS PHÉNOLOGIQUES
# ==============================================================================
if (!SKIP_VIZ) {
  cli::cli_h2("Étape 2 — Visualisation des profils phénologiques")

  p_profiles <- plot_phenological_profiles(
    ts_long,
    index_name = "NDVI",
    save_path  = file.path(FIGURES_DIR, "01_phenological_profiles_NDVI.png")
  )

  plot_deciduous_vs_evergreen(
    ts_long,
    save_path = file.path(FIGURES_DIR, "02_deciduous_vs_evergreen.png")
  )

  plot_profiles_by_type(
    ts_long,
    save_path = file.path(FIGURES_DIR, "03_profiles_by_type.png")
  )

  cli::cli_alert_success("Graphiques phénologiques générés")
}

# ==============================================================================
# ÉTAPE 3 : CONSTRUCTION DE LA MATRICE DE FEATURES
# ==============================================================================
cli::cli_h2("Étape 3 — Construction de la matrice de features")

if (use_multiyear && USE_SYNTHETIC) {
  cli::cli_alert_info("Construction multi-année ({length(YEARS)} ans)")
  feature_matrix <- build_multiyear_features(
    years                 = YEARS,
    n_samples_per_species = N_SAMPLES
  )
} else {
  feature_matrix <- build_feature_matrix(ts_long)
}

fm_path <- file.path(PROCESSED_DIR, "feature_matrix.csv")
readr::write_csv(feature_matrix, fm_path)
cli::cli_alert_success("Matrice de features sauvegardée : {fm_path}")
cli::cli_text("  Dimensions : {nrow(feature_matrix)} × {ncol(feature_matrix)}")

# ==============================================================================
# ÉTAPE 4 : SÉPARATION TRAIN / TEST
# ==============================================================================
cli::cli_h2("Étape 4 — Séparation train/test")

split <- split_train_test(feature_matrix, target_col = "species_name")
train_data <- split$train
test_data  <- split$test

feature_cols <- select_features(feature_matrix)
cli::cli_text("  {length(feature_cols)} features sélectionnées (règles)")

# ==============================================================================
# ÉTAPE 4b : SÉLECTION DE FEATURES PAR BORUTA (optionnel)
# ==============================================================================
boruta_result <- NULL
if (isTRUE(USE_BORUTA)) {
  cli::cli_h2("Étape 4b — Sélection de features (Boruta)")

  boruta_out <- select_features_boruta(
    feature_matrix = train_data,
    feature_cols   = feature_cols,
    max_runs       = BORUTA_MAXRUNS
  )

  boruta_result <- boruta_out$boruta_result
  feature_cols  <- boruta_out$selected_cols

  cli::cli_text("  {length(feature_cols)} features retenues après Boruta")

  if (!is.null(boruta_out$importance_df)) {
    readr::write_csv(boruta_out$importance_df,
                     file.path(OUTPUT_DIR, "boruta_importance.csv"))
  }
}

# ==============================================================================
# ÉTAPE 5 : CLASSIFICATION
# ==============================================================================
cli::cli_h2("Étape 5 — Classification")

if (MODE == "rf" || MODE == "both") {
  cli::cli_h3("5a. Random Forest (ranger)")

  rf_model <- train_random_forest(train_data, feature_cols)
  rf_preds <- predict_rf(rf_model, test_data)

  rf_eval <- evaluate_classification(
    y_true      = test_data$species_name,
    y_pred      = rf_preds$predicted_class,
    class_names = eval_class_names
  )

  rf_importance <- get_variable_importance(rf_model)
  save_model(rf_model, rf_eval, model_name = "treesatai_rf")

  readr::write_csv(rf_importance, file.path(OUTPUT_DIR, "rf_variable_importance.csv"))
  readr::write_csv(rf_eval$per_class, file.path(OUTPUT_DIR, "rf_per_class_metrics.csv"))
}

if (MODE == "cnn" || MODE == "both") {
  cli::cli_h3("5b. CNN temporel (torch)")

  if (requireNamespace("torch", quietly = TRUE)) {
    train_tensors <- prepare_cnn_data(train_data)
    test_tensors  <- prepare_cnn_data(test_data)

    if (!is.null(train_tensors)) {
      cnn_result <- train_temporal_cnn(
        train_tensors,
        val_tensors = test_tensors,
        n_epochs    = 50,
        batch_size  = 32,
        lr          = 0.001
      )

      cnn_result$model$eval()
      torch::with_no_grad({
        test_output <- cnn_result$model(test_tensors$X)
        cnn_preds_idx <- test_output$argmax(dim = 2)$to(device = "cpu")$numpy()
      })
      cnn_preds <- cnn_result$class_names[cnn_preds_idx + 1]

      cnn_eval <- evaluate_classification(
        y_true      = test_data$species_name,
        y_pred      = cnn_preds,
        class_names = eval_class_names
      )

      readr::write_csv(cnn_eval$per_class, file.path(OUTPUT_DIR, "cnn_per_class_metrics.csv"))
      readr::write_csv(cnn_result$history, file.path(OUTPUT_DIR, "cnn_training_history.csv"))
    }
  } else {
    cli::cli_alert_warning("torch non installé — CNN ignoré")
  }
}

# ==============================================================================
# ÉTAPE 6 : VALIDATION CROISÉE
# ==============================================================================
cli::cli_h2("Étape 6 — Validation croisée")

cv_results <- cross_validate(feature_matrix, feature_cols,
                              k = CLASSIF_PARAMS$cv_folds,
                              repeats = CLASSIF_PARAMS$cv_repeats)

readr::write_csv(cv_results$results, file.path(OUTPUT_DIR, "cv_results.csv"))

# ==============================================================================
# ÉTAPE 7 : VISUALISATIONS FINALES
# ==============================================================================
if (!SKIP_VIZ) {
  cli::cli_h2("Étape 7 — Visualisations finales")

  eval_results <- if (exists("rf_eval")) rf_eval else cnn_eval

  plot_confusion_matrix(
    eval_results$confusion_matrix,
    class_names = eval_class_names,
    title = glue::glue("Matrice de confusion — {n_classes} classes"),
    save_path = file.path(FIGURES_DIR, "04_confusion_matrix.png")
  )

  if (exists("rf_importance")) {
    plot_variable_importance(
      rf_importance,
      top_n = 30,
      save_path = file.path(FIGURES_DIR, "05_variable_importance.png")
    )
  }

  plot_species_metrics(
    eval_results$per_class,
    save_path = file.path(FIGURES_DIR, "06_species_metrics.png")
  )

  plot_spectral_heatmap(
    feature_matrix,
    save_path = file.path(FIGURES_DIR, "07_spectral_heatmap.png")
  )

  if (exists("rf_importance")) {
    dashboard_plots <- list(
      profiles   = if (exists("p_profiles")) p_profiles else NULL,
      confusion  = NULL,
      importance = NULL,
      metrics    = NULL
    )
    plot_dashboard(
      dashboard_plots,
      save_path = file.path(FIGURES_DIR, "08_dashboard.png")
    )
  }

  cli::cli_alert_success("Toutes les visualisations générées dans {FIGURES_DIR}")
}

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================
cli::cli_h1("Résumé de l'exécution")
cli::cli_text("")

cli::cli_alert_success("Pipeline terminé avec succès !")
cli::cli_text("")

summary_table <- data.frame(
  Paramètre = c(
    "Classes", "Échantillons", "Features",
    "Méthode", "Multi-année", "Boruta",
    "OA (%)", "Kappa", "Macro F1 (%)"
  ),
  Valeur = c(
    n_classes,
    nrow(feature_matrix),
    length(feature_cols),
    toupper(MODE),
    if (use_multiyear) paste(YEARS, collapse=", ") else "Non",
    if (!is.null(boruta_result)) "Oui" else "Non",
    if (exists("rf_eval")) round(rf_eval$overall_accuracy * 100, 2) else
      if (exists("cnn_eval")) round(cnn_eval$overall_accuracy * 100, 2) else "N/A",
    if (exists("rf_eval")) round(rf_eval$kappa, 4) else
      if (exists("cnn_eval")) round(cnn_eval$kappa, 4) else "N/A",
    if (exists("rf_eval")) round(rf_eval$macro_f1 * 100, 2) else
      if (exists("cnn_eval")) round(cnn_eval$macro_f1 * 100, 2) else "N/A"
  )
)

print(summary_table, row.names = FALSE)

cli::cli_text("")
cli::cli_text("Fichiers de sortie :")
cli::cli_ul()
cli::cli_li("Données     : {PROCESSED_DIR}")
cli::cli_li("Modèles     : {MODELS_DIR}")
cli::cli_li("Figures     : {FIGURES_DIR}")
cli::cli_li("Résultats   : {OUTPUT_DIR}")
cli::cli_end()
