# ==============================================================================
# TreeSatAI-Time-Series (IGNF) — Pipeline complet (fonction package)
# Classification de 10 classes forestières européennes
# (20 espèces regroupées) par séries temporelles Sentinel-1/2 annuelles
#
# Usage depuis R :
#   train_treesatai()                              # Données synthétiques
#   train_treesatai(data_path = "/chemin/dataset") # Données réelles
#   train_treesatai(n_samples = 100)               # Plus d'échantillons
#
# Usage CLI (via inst/scripts/06_pipeline.R) :
#   Rscript inst/scripts/06_pipeline.R --synthetic
#   Rscript inst/scripts/06_pipeline.R --data /chemin/dataset
# ==============================================================================

#' Entraîner un modèle TreeSatAI sur données synthétiques ou réelles
#'
#' @param data_path Chemin vers un répertoire de données réelles (NULL = synthétique)
#' @param mode Méthode : "rf" (Random Forest), "cnn" (CNN temporel), ou "both"
#' @param n_samples Nombre d'échantillons par espèce en mode synthétique
#' @param year Année de simulation (mode synthétique)
#' @param skip_viz Ne pas générer les visualisations
#' @param output_dir Répertoire de sortie
#' @return Liste avec le modèle, l'évaluation et la matrice de features
#' @export
train_treesatai <- function(data_path   = NULL,
                             mode        = "rf",
                             n_samples   = 50,
                             year        = 2021,
                             skip_viz    = FALSE,
                             output_dir  = NULL) {

  root <- .get_project_root()
  if (is.null(output_dir)) output_dir <- file.path(root, "output")
  figures_dir   <- file.path(root, "figures")
  processed_dir <- file.path(root, "data", "processed")
  models_dir    <- file.path(output_dir, "models")

  init_project_dirs(root)

  # ===========================================================================
  cli::cli_h1("Pipeline TreeSatAI-Time-Series")
  cli::cli_text("")
  # ===========================================================================

  # ===========================================================================
  # ÉTAPE 1 : ACQUISITION DES DONNÉES
  # ===========================================================================
  cli::cli_h2("\u00c9tape 1 \u2014 Acquisition des donn\u00e9es")

  if (is.null(data_path)) {
    # ---- Mode données synthétiques TreeSatAI ----
    cli::cli_alert_info("Mode synth\u00e9tique : g\u00e9n\u00e9ration de {n_samples} \u00e9chantillons \u00d7 {nrow(SPECIES)} esp\u00e8ces")
    cli::cli_text("")
    cli::cli_text("Les profils ph\u00e9nologiques sont bas\u00e9s sur les signatures spectrales")
    cli::cli_text("caract\u00e9ristiques de chaque essence (double logistique + bruit).")
    cli::cli_text("")

    ts_long <- generate_synthetic_dataset(n_samples_per_species = n_samples, year = year)
    cli::cli_alert_success("Dataset synth\u00e9tique g\u00e9n\u00e9r\u00e9 : {nrow(ts_long)} observations")
    use_treesatai_split <- FALSE
    treesatai <- NULL

  } else {
    # ---- Mode données réelles ----
    cli::cli_alert_info("Chargement des donn\u00e9es depuis : {data_path}")

    if (!dir.exists(data_path)) {
      cli::cli_alert_danger("R\u00e9pertoire introuvable : {data_path}")
      stop("Chemin invalide")
    }

    # Détecter le format TreeSatAI (labels JSON multi-genres)
    labels_file <- list.files(data_path, pattern = "multi_labels.*\\.json$",
                               recursive = TRUE, full.names = TRUE)

    if (length(labels_file) > 0) {
      # ---- Mode TreeSatAI (patches + labels JSON) ----
      treesatai <- load_treesatai_data(data_path)
      ts_long <- treesatai$ts_long
      treesatai_split <- treesatai$split
      use_treesatai_split <- !is.null(treesatai_split)

    } else {
      # ---- Mode générique (plots shapefile + Sentinel-2 scenes) ----
      use_treesatai_split <- FALSE
      treesatai <- NULL

      plot_files <- list.files(data_path, pattern = "\\.(gpkg|shp|geojson)$",
                               recursive = TRUE, full.names = TRUE)
      plots <- NULL
      if (length(plot_files) > 0) {
        plots <- sf::st_read(plot_files[1], quiet = TRUE)
        cli::cli_alert_success("{nrow(plots)} parcelles charg\u00e9es")
      }

      ts_files <- list.files(data_path, pattern = "\\.(csv|parquet)$",
                             recursive = TRUE, full.names = TRUE)
      if (length(ts_files) > 0) {
        ts_long <- readr::read_csv(ts_files[1], show_col_types = FALSE)
        cli::cli_alert_success("S\u00e9ries temporelles charg\u00e9es : {nrow(ts_long)} lignes")
      } else if (!is.null(plots)) {
        s2_dir <- list.dirs(data_path, recursive = FALSE)
        s2_dir <- s2_dir[grep("sentinel|S2", s2_dir, ignore.case = TRUE)]
        if (length(s2_dir) > 0) {
          ts_long <- extract_s2_timeseries(plots, s2_dir[1])
        } else {
          cli::cli_alert_danger("Aucune donn\u00e9e Sentinel-2 trouv\u00e9e dans {data_path}")
          stop("Donn\u00e9es manquantes")
        }
      } else {
        cli::cli_alert_danger("Aucune parcelle ni s\u00e9rie temporelle trouv\u00e9e dans {data_path}")
        stop("Donn\u00e9es manquantes : placez des .gpkg/.shp/.geojson ou des .csv/.parquet")
      }
    }
  }

  # Regroupement des 20 espèces en 10 classes forestières
  cli::cli_alert_info("Regroupement des esp\u00e8ces en {length(unique(SPECIES_GROUPS))} classes")
  ts_long <- remap_species_groups(ts_long)

  # Déterminer les noms de classes pour l'évaluation
  eval_class_names <- if (!is.null(treesatai) && !is.null(treesatai$genus_names)) {
    # Remapper aussi les noms TreeSatAI si applicable
    gnames <- treesatai$genus_names
    gnames <- ifelse(gnames %in% names(SPECIES_GROUPS), SPECIES_GROUPS[gnames], gnames)
    sort(unique(gnames))
  } else {
    SPECIES_GROUPS_INFO$group
  }

  # ===========================================================================
  # ÉTAPE 2 : VISUALISATION DES PROFILS PHÉNOLOGIQUES
  # ===========================================================================
  p_profiles <- NULL
  if (!skip_viz) {
    cli::cli_h2("\u00c9tape 2 \u2014 Visualisation des profils ph\u00e9nologiques")

    p_profiles <- plot_phenological_profiles(
      ts_long,
      index_name = "NDVI",
      save_path  = file.path(figures_dir, "01_phenological_profiles_NDVI.png")
    )

    plot_deciduous_vs_evergreen(
      ts_long,
      save_path = file.path(figures_dir, "02_deciduous_vs_evergreen.png")
    )

    plot_profiles_by_type(
      ts_long,
      save_path = file.path(figures_dir, "03_profiles_by_type.png")
    )

    cli::cli_alert_success("Graphiques ph\u00e9nologiques g\u00e9n\u00e9r\u00e9s")
  }

  # ===========================================================================
  # ÉTAPE 3 : CONSTRUCTION DE LA MATRICE DE FEATURES
  # ===========================================================================
  cli::cli_h2("\u00c9tape 3 \u2014 Construction de la matrice de features")

  feature_matrix <- build_feature_matrix(ts_long)

  fm_path <- file.path(processed_dir, "feature_matrix.csv")
  readr::write_csv(feature_matrix, fm_path)
  cli::cli_alert_success("Matrice de features sauvegard\u00e9e : {fm_path}")
  cli::cli_text("  Dimensions : {nrow(feature_matrix)} \u00d7 {ncol(feature_matrix)}")

  # ===========================================================================
  # ÉTAPE 4 : SÉPARATION TRAIN / TEST
  # ===========================================================================
  cli::cli_h2("\u00c9tape 4 \u2014 S\u00e9paration train/test")

  if (exists("use_treesatai_split") && isTRUE(use_treesatai_split)) {
    # Utiliser le split prédéfini du dataset TreeSatAI
    cli::cli_alert_info("Utilisation du split pr\u00e9d\u00e9fini TreeSatAI")
    train_patches <- treesatai_split$train
    test_patches  <- treesatai_split$test

    train_data <- feature_matrix[feature_matrix$plot_id %in% train_patches, ]
    test_data  <- feature_matrix[feature_matrix$plot_id %in% test_patches, ]

    # Patchs non trouvés dans le split → ajouter au train
    unmatched <- feature_matrix[!feature_matrix$plot_id %in% c(train_patches, test_patches), ]
    if (nrow(unmatched) > 0) {
      cli::cli_alert_warning("{nrow(unmatched)} patchs non trouv\u00e9s dans le split, ajout\u00e9s au train")
      train_data <- rbind(train_data, unmatched)
    }

    cli::cli_alert_success("Train : {nrow(train_data)} / Test : {nrow(test_data)}")
  } else {
    split <- split_train_test(feature_matrix, target_col = "species_name")
    train_data <- split$train
    test_data  <- split$test
  }

  feature_cols <- select_features(feature_matrix)
  cli::cli_text("  {length(feature_cols)} features s\u00e9lectionn\u00e9es")

  # ===========================================================================
  # ÉTAPE 5 : CLASSIFICATION
  # ===========================================================================
  cli::cli_h2("\u00c9tape 5 \u2014 Classification")

  rf_model <- NULL
  rf_eval  <- NULL
  rf_importance <- NULL
  cnn_eval <- NULL

  if (mode == "rf" || mode == "both") {
    cli::cli_h3("5a. Random Forest (ranger)")

    rf_model <- train_random_forest(train_data, feature_cols)
    rf_preds <- predict_rf(rf_model, test_data)

    rf_eval <- evaluate_classification(
      y_true      = test_data$species_name,
      y_pred      = rf_preds$predicted_class,
      class_names = eval_class_names
    )

    rf_importance <- get_variable_importance(rf_model)
    save_model(rf_model, rf_eval, model_name = "treesatai_rf", output_dir = models_dir)

    readr::write_csv(rf_importance, file.path(output_dir, "rf_variable_importance.csv"))
    readr::write_csv(rf_eval$per_class, file.path(output_dir, "rf_per_class_metrics.csv"))
  }

  if (mode == "cnn" || mode == "both") {
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

        readr::write_csv(cnn_eval$per_class, file.path(output_dir, "cnn_per_class_metrics.csv"))
        readr::write_csv(cnn_result$history, file.path(output_dir, "cnn_training_history.csv"))
      }
    } else {
      cli::cli_alert_warning("torch non install\u00e9 \u2014 CNN ignor\u00e9")
    }
  }

  # ===========================================================================
  # ÉTAPE 6 : VALIDATION CROISÉE
  # ===========================================================================
  cli::cli_h2("\u00c9tape 6 \u2014 Validation crois\u00e9e")

  cv_results <- cross_validate(feature_matrix, feature_cols,
                                k = CLASSIF_PARAMS$cv_folds,
                                repeats = CLASSIF_PARAMS$cv_repeats)

  readr::write_csv(cv_results$results, file.path(output_dir, "cv_results.csv"))

  # ===========================================================================
  # ÉTAPE 7 : VISUALISATIONS FINALES
  # ===========================================================================
  if (!skip_viz) {
    cli::cli_h2("\u00c9tape 7 \u2014 Visualisations finales")

    eval_results <- if (!is.null(rf_eval)) rf_eval else cnn_eval

    n_classes <- length(eval_class_names)
    plot_confusion_matrix(
      eval_results$confusion_matrix,
      class_names = eval_class_names,
      title = glue::glue("Matrice de confusion \u2014 {n_classes} classes"),
      save_path = file.path(figures_dir, "04_confusion_matrix.png")
    )

    if (!is.null(rf_importance)) {
      plot_variable_importance(
        rf_importance,
        top_n = 30,
        save_path = file.path(figures_dir, "05_variable_importance.png")
      )
    }

    plot_species_metrics(
      eval_results$per_class,
      save_path = file.path(figures_dir, "06_species_metrics.png")
    )

    plot_spectral_heatmap(
      feature_matrix,
      save_path = file.path(figures_dir, "07_spectral_heatmap.png")
    )

    cli::cli_alert_success("Toutes les visualisations g\u00e9n\u00e9r\u00e9es dans {figures_dir}")
  }

  # ===========================================================================
  # RÉSUMÉ FINAL
  # ===========================================================================
  eval_final <- if (!is.null(rf_eval)) rf_eval else cnn_eval

  cli::cli_h1("R\u00e9sum\u00e9")
  cli::cli_alert_success("Pipeline termin\u00e9 avec succ\u00e8s !")
  if (!is.null(eval_final)) {
    cli::cli_alert_info("OA : {round(eval_final$overall_accuracy * 100, 1)}%  |  Kappa : {round(eval_final$kappa, 3)}")
  }
  cli::cli_text("")
  cli::cli_text("Mod\u00e8le sauvegard\u00e9 dans : {.path {models_dir}}")

  invisible(list(
    model          = rf_model,
    evaluation     = rf_eval,
    feature_matrix = feature_matrix,
    feature_cols   = feature_cols,
    cv_results     = cv_results
  ))
}
