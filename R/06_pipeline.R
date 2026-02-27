#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series (IGNF) — Pipeline complet
# Classification de 20 essences forestières européennes
# par séries temporelles Sentinel-2 annuelles
#
# Usage :
#   Rscript R/06_pipeline.R --data /chemin/dataset  # Données réelles
#   Rscript R/06_pipeline.R --data /chemin --mode cnn  # Avec CNN temporel
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

DATA_PATH   <- NULL
MODE        <- "rf"       # "rf" ou "cnn"
SKIP_VIZ    <- FALSE

for (i in seq_along(args)) {
  if (args[i] == "--data" && i < length(args))    DATA_PATH <- args[i + 1]
  if (args[i] == "--mode" && i < length(args))    MODE      <- args[i + 1]
  if (args[i] == "--no-viz")                       SKIP_VIZ  <- TRUE
}

# ==============================================================================
cli::cli_h1("Pipeline TreeSatAI-Time-Series")
cli::cli_text("")
# ==============================================================================

# ==============================================================================
# ÉTAPE 1 : ACQUISITION DES DONNÉES
# ==============================================================================
cli::cli_h2("Étape 1 — Acquisition des données")

if (!is.null(DATA_PATH)) {
  # ---- Mode données réelles ----
  cli::cli_alert_info("Chargement des données depuis : {DATA_PATH}")

  if (dir.exists(DATA_PATH)) {
    # Rechercher le fichier de parcelles
    plot_files <- list.files(DATA_PATH, pattern = "\\.(gpkg|shp|geojson)$",
                             recursive = TRUE, full.names = TRUE)
    if (length(plot_files) > 0) {
      plots <- sf::st_read(plot_files[1], quiet = TRUE)
      cli::cli_alert_success("{nrow(plots)} parcelles chargées")
    }

    # Rechercher les séries temporelles pré-extraites
    ts_files <- list.files(DATA_PATH, pattern = "\\.(csv|parquet)$",
                           recursive = TRUE, full.names = TRUE)
    if (length(ts_files) > 0) {
      ts_long <- readr::read_csv(ts_files[1], show_col_types = FALSE)
      cli::cli_alert_success("Séries temporelles chargées : {nrow(ts_long)} lignes")
    } else {
      # Extraire les séries temporelles depuis les rasters Sentinel-2
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
  cli::cli_text("  Rscript R/06_pipeline.R --data /chemin/vers/donnees")
  cli::cli_text("")
  cli::cli_text("Le répertoire doit contenir :")
  cli::cli_ul()
  cli::cli_li("Un fichier vectoriel (.gpkg, .shp) avec les parcelles et espèces")
  cli::cli_li("Des séries temporelles (.csv) ou des images Sentinel-2 brutes")
  cli::cli_end()
  stop("Argument --data obligatoire. Pas de données, pas de pipeline.")
}

# ==============================================================================
# ÉTAPE 2 : VISUALISATION DES PROFILS PHÉNOLOGIQUES
# ==============================================================================
if (!SKIP_VIZ) {
  cli::cli_h2("Étape 2 — Visualisation des profils phénologiques")

  # Profils NDVI de toutes les espèces
  p_profiles <- plot_phenological_profiles(
    ts_long,
    index_name = "NDVI",
    save_path  = file.path(FIGURES_DIR, "01_phenological_profiles_NDVI.png")
  )

  # Comparaison caduc vs persistant
  p_dec_ever <- plot_deciduous_vs_evergreen(
    ts_long,
    save_path = file.path(FIGURES_DIR, "02_deciduous_vs_evergreen.png")
  )

  # Profils par type
  p_types <- plot_profiles_by_type(
    ts_long,
    save_path = file.path(FIGURES_DIR, "03_profiles_by_type.png")
  )

  cli::cli_alert_success("Graphiques phénologiques générés")
}

# ==============================================================================
# ÉTAPE 3 : CONSTRUCTION DE LA MATRICE DE FEATURES
# ==============================================================================
cli::cli_h2("Étape 3 — Construction de la matrice de features")

feature_matrix <- build_feature_matrix(ts_long)

# Sauvegarder la matrice
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

# Sélection des features
feature_cols <- select_features(feature_matrix)
cli::cli_text("  {length(feature_cols)} features sélectionnées")

# ==============================================================================
# ÉTAPE 5 : CLASSIFICATION
# ==============================================================================
cli::cli_h2("Étape 5 — Classification")

if (MODE == "rf" || MODE == "both") {
  # ---- Random Forest ----
  cli::cli_h3("5a. Random Forest (ranger)")

  rf_model <- train_random_forest(train_data, feature_cols)

  # Prédiction
  rf_preds <- predict_rf(rf_model, test_data)

  # Évaluation
  rf_eval <- evaluate_classification(
    y_true      = test_data$species_name,
    y_pred      = rf_preds$predicted_class,
    class_names = SPECIES$french
  )

  # Importance des variables
  rf_importance <- get_variable_importance(rf_model)

  # Sauvegarder
  save_model(rf_model, rf_eval, model_name = "treesatai_rf")

  # Export des résultats
  readr::write_csv(rf_importance, file.path(OUTPUT_DIR, "rf_variable_importance.csv"))
  readr::write_csv(rf_eval$per_class, file.path(OUTPUT_DIR, "rf_per_class_metrics.csv"))
}

if (MODE == "cnn" || MODE == "both") {
  # ---- CNN temporel ----
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

      # Prédiction finale
      cnn_result$model$eval()
      torch::with_no_grad({
        test_output <- cnn_result$model(test_tensors$X)
        cnn_preds_idx <- test_output$argmax(dim = 2)$to(device = "cpu")$numpy()
      })
      cnn_preds <- cnn_result$class_names[cnn_preds_idx + 1]

      # Évaluation
      cnn_eval <- evaluate_classification(
        y_true      = test_data$species_name,
        y_pred      = cnn_preds,
        class_names = SPECIES$french
      )

      # Sauvegarder
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

  # Utiliser les résultats RF par défaut
  eval_results <- if (exists("rf_eval")) rf_eval else cnn_eval

  # Matrice de confusion
  p_conf <- plot_confusion_matrix(
    eval_results$confusion_matrix,
    class_names = SPECIES$french,
    title = "Matrice de confusion — 20 espèces",
    save_path = file.path(FIGURES_DIR, "04_confusion_matrix.png")
  )

  # Importance des variables (RF uniquement)
  if (exists("rf_importance")) {
    p_imp <- plot_variable_importance(
      rf_importance,
      top_n = 30,
      save_path = file.path(FIGURES_DIR, "05_variable_importance.png")
    )
  }

  # Métriques par espèce
  p_metrics <- plot_species_metrics(
    eval_results$per_class,
    save_path = file.path(FIGURES_DIR, "06_species_metrics.png")
  )

  # Heatmap des signatures
  plot_spectral_heatmap(
    feature_matrix,
    save_path = file.path(FIGURES_DIR, "07_spectral_heatmap.png")
  )

  # Dashboard récapitulatif
  if (exists("rf_importance") && exists("p_imp")) {
    dashboard_plots <- list(
      profiles   = p_profiles,
      confusion  = p_conf,
      importance = p_imp,
      metrics    = p_metrics
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

# Tableau récapitulatif
summary_table <- data.frame(
  Paramètre = c(
    "Espèces", "Échantillons", "Features",
    "Méthode", "OA (%)", "Kappa", "Macro F1 (%)"
  ),
  Valeur = c(
    nrow(SPECIES),
    nrow(feature_matrix),
    length(feature_cols),
    toupper(MODE),
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

cli::cli_text("")
cli::cli_text("Profils phénologiques caractéristiques :")
cli::cli_ul()
cli::cli_li("Feuillus caducifoliés : forte amplitude NDVI (0.15–0.85), pic en juin-juillet")
cli::cli_li("Résineux persistants  : NDVI stable (0.45–0.60), faible variation saisonnière")
cli::cli_li("Mélèze (résineux caduc) : profil intermédiaire, perte d'aiguilles en automne")
cli::cli_li("Chêne vert (feuillu persistant) : NDVI stable et élevé toute l'année")
cli::cli_end()
