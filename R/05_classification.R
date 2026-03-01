#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Classification d'essences forestières
# Random Forest + CNN temporel optionnel
# ==============================================================================

# Configuration, utilitaires et phénologie chargés via le package

# ==============================================================================
# 1. PRÉPARATION DES DONNÉES POUR LA CLASSIFICATION
# ==============================================================================

#' Séparation train/test stratifiée par espèce
#' @param feature_matrix data.frame de features (1 ligne par parcelle)
#' @param target_col Nom de la colonne cible (espèce)
#' @param train_ratio Proportion d'entraînement
#' @param seed Graine aléatoire
#' @return Liste avec train et test
split_train_test <- function(feature_matrix, target_col = "species_name",
                              train_ratio = CLASSIF_PARAMS$train_ratio,
                              seed = CLASSIF_PARAMS$seed) {
  set.seed(seed)
  log_msg("Séparation train/test ({train_ratio*100}% / {(1-train_ratio)*100}%)")

  species <- unique(feature_matrix[[target_col]])
  train_idx <- c()

  for (sp in species) {
    sp_idx <- which(feature_matrix[[target_col]] == sp)
    n_train <- round(length(sp_idx) * train_ratio)
    selected <- sample(sp_idx, n_train)
    train_idx <- c(train_idx, selected)
  }

  train_data <- feature_matrix[train_idx, ]
  test_data  <- feature_matrix[-train_idx, ]

  log_msg("  Train : {nrow(train_data)} parcelles", level = "info")
  log_msg("  Test  : {nrow(test_data)} parcelles", level = "info")

  # Vérifier l'équilibre des classes
  train_counts <- table(train_data[[target_col]])
  test_counts  <- table(test_data[[target_col]])

  log_msg("  Distribution train : min={min(train_counts)}, max={max(train_counts)}, median={median(train_counts)}",
          level = "info")

  list(train = train_data, test = test_data)
}

#' Sélection des features à utiliser pour la classification
#' @param feature_matrix data.frame de features
#' @return Vecteur des noms de colonnes features
select_features <- function(feature_matrix) {
  all_cols <- names(feature_matrix)

  # Colonnes à exclure (identifiants)
  exclude <- c("plot_id", "species_code", "species_name")

  # Sélection selon la configuration
  feature_cols <- c()

  if (CLASSIF_PARAMS$use_spectral_bands) {
    # Bandes temporelles : B02_d001, B02_d006, ...
    band_cols <- grep("^B\\d{2}_d\\d+|^B8A_d\\d+", all_cols, value = TRUE)
    feature_cols <- c(feature_cols, band_cols)
  }

  if (CLASSIF_PARAMS$use_spectral_indices) {
    # Indices temporels : NDVI_d001, EVI_d001, ...
    index_cols <- grep("^(NDVI|EVI|NDWI|CRI|RENDVI|NBR)_d\\d+", all_cols, value = TRUE)
    # Statistiques des indices : NDVI_mean, NDVI_sd, ...
    stats_cols <- grep("^(NDVI|EVI|NDWI|CRI|RENDVI|NBR)_(mean|sd|min|max|range|median|q25|q75|iqr)",
                       all_cols, value = TRUE)
    feature_cols <- c(feature_cols, index_cols, stats_cols)
  }

  if (CLASSIF_PARAMS$use_phenometrics) {
    # Métriques phénologiques : pheno_SOS, pheno_EOS, etc.
    pheno_cols <- grep("pheno_", all_cols, value = TRUE)
    feature_cols <- c(feature_cols, pheno_cols)
  }

  if (CLASSIF_PARAMS$use_temporal_stats) {
    # Statistiques des bandes brutes : B02_mean, B02_sd, ...
    band_stats <- grep("^B\\d{2}_(mean|sd|min|max|range|median|q25|q75|iqr)", all_cols, value = TRUE)
    b8a_stats  <- grep("^B8A_(mean|sd|min|max|range|median|q25|q75|iqr)", all_cols, value = TRUE)
    feature_cols <- c(feature_cols, band_stats, b8a_stats)
  }

  # Fourier
  fourier_cols <- grep("fourier_", all_cols, value = TRUE)
  feature_cols <- c(feature_cols, fourier_cols)

  # Features topographiques (MNT, pente, exposition, TWI, TPI)
  if (isTRUE(CLASSIF_PARAMS$use_terrain)) {
    terrain_cols <- grep("^DEM_", all_cols, value = TRUE)
    feature_cols <- c(feature_cols, terrain_cols)
  }

  # Déduplication
  feature_cols <- unique(feature_cols)

  # Retirer les colonnes exclues
  feature_cols <- setdiff(feature_cols, exclude)

  # Retirer les colonnes avec trop de NA
  na_ratio <- sapply(feature_matrix[feature_cols], function(x) mean(is.na(x)))
  feature_cols <- feature_cols[na_ratio < 0.5]

  log_msg("Features sélectionnées : {length(feature_cols)} colonnes")

  feature_cols
}

#' Sélection de features par Boruta (wrapper Random Forest)
#'
#' Boruta crée des « shadow features » (copies aléatoirement permutées)
#' et teste si les features réelles surpassent les shadows. Résultat :
#' Confirmed / Rejected / Tentative.
#'
#' @param feature_matrix data.frame avec colonnes features + species_name
#' @param feature_cols Vecteur des noms de features candidates
#' @param target_col Colonne cible
#' @param max_runs Nombre max d'itérations Boruta (défaut : 100)
#' @param p_value Seuil de significativité (défaut : 0.01)
#' @return Liste avec selected_cols, boruta_result, importance_df
#' @export
select_features_boruta <- function(feature_matrix,
                                    feature_cols,
                                    target_col = "species_name",
                                    max_runs   = 100,
                                    p_value    = 0.01) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    cli::cli_alert_warning("{.pkg Boruta} non install\u00e9 \u2014 s\u00e9lection Boruta ignor\u00e9e")
    cli::cli_text("Installer avec : {.code install.packages('Boruta')}")
    return(list(selected_cols = feature_cols, boruta_result = NULL, importance_df = NULL))
  }

  log_msg("S\u00e9lection de features par Boruta ({length(feature_cols)} candidates)...")
  log_msg("  maxRuns={max_runs}, pValue={p_value}", level = "info")

  # Préparer les données
  X <- feature_matrix[, feature_cols, drop = FALSE]
  y <- as.factor(feature_matrix[[target_col]])

  # Remplacement NA par médiane
  for (col in feature_cols) {
    na_mask <- is.na(X[[col]])
    if (any(na_mask)) {
      X[[col]][na_mask] <- median(X[[col]], na.rm = TRUE)
    }
  }

  # Retirer colonnes à variance nulle (Boruta échoue sinon)
  var_check <- sapply(X, var, na.rm = TRUE)
  zero_var <- names(var_check[var_check == 0 | is.na(var_check)])
  if (length(zero_var) > 0) {
    X <- X[, !names(X) %in% zero_var, drop = FALSE]
    log_msg("  {length(zero_var)} features \u00e0 variance nulle retir\u00e9es", level = "warning")
  }

  # Exécuter Boruta avec ranger comme backend (rapide)
  boruta_result <- Boruta::Boruta(
    x        = X,
    y        = y,
    pValue   = p_value,
    maxRuns  = max_runs,
    doTrace  = 0,
    num.trees = 200,
    num.threads = max(1, parallel::detectCores() - 1)
  )

  # Résoudre les features tentatives (conservateur : les garder)
  boruta_final <- Boruta::TentativeRoughFix(boruta_result)

  # Extraire les features confirmées + tentatives résolues
  decision <- boruta_final$finalDecision
  confirmed <- names(decision[decision == "Confirmed"])
  rejected  <- names(decision[decision == "Rejected"])
  tentative <- names(decision[decision == "Tentative"])

  log_msg("Boruta termin\u00e9 :", level = "success")
  log_msg("  Confirm\u00e9es  : {length(confirmed)}", level = "info")
  log_msg("  Rejet\u00e9es   : {length(rejected)}", level = "info")
  if (length(tentative) > 0) {
    log_msg("  Tentatives : {length(tentative)} (conserv\u00e9es)", level = "info")
  }

  selected_cols <- c(confirmed, tentative)

  # Tableau d'importance pour diagnostic
  imp_df <- data.frame(
    variable   = names(decision),
    decision   = as.character(decision),
    meanImp    = apply(boruta_result$ImpHistory, 2, mean, na.rm = TRUE)[names(decision)],
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(-imp_df$meanImp), ]

  log_msg("  R\u00e9duction : {length(feature_cols)} \u2192 {length(selected_cols)} features ({round(length(selected_cols)/length(feature_cols)*100)}%)")

  list(
    selected_cols  = selected_cols,
    boruta_result  = boruta_final,
    importance_df  = imp_df
  )
}

# ==============================================================================
# 2. RANDOM FOREST
# ==============================================================================

#' Calcul des poids de classe pour compenser le déséquilibre
#'
#' Retourne un vecteur de poids par échantillon (case weights) inversement
#' proportionnels à la fréquence de chaque classe : w_i = N / (K × n_i)
#' où N = effectif total, K = nombre de classes, n_i = effectif de la classe i.
#'
#' @param y Vecteur factor ou character des labels de classe
#' @return Vecteur numérique de même longueur que y (un poids par échantillon)
#' @export
compute_class_weights <- function(y) {
  y <- as.factor(y)
  class_counts <- table(y)
  n_total   <- length(y)
  n_classes <- length(class_counts)

  # Poids par classe : inversement proportionnel à la fréquence
  weight_per_class <- n_total / (n_classes * class_counts)

  # Attribuer à chaque échantillon le poids de sa classe
  sample_weights <- as.numeric(weight_per_class[as.character(y)])

  log_msg("Poids de classe (déséquilibre) :")
  for (cls in names(sort(class_counts))) {
    log_msg("  {cls} : n={class_counts[cls]}, poids={round(weight_per_class[cls], 3)}",
            level = "info")
  }

  sample_weights
}

#' Entraînement d'un Random Forest pour la classification d'espèces
#' @param train_data data.frame d'entraînement
#' @param feature_cols Colonnes features à utiliser
#' @param target_col Colonne cible
#' @return Modèle randomForest ajusté
train_random_forest <- function(train_data, feature_cols, target_col = "species_name") {
  log_msg("Entraînement Random Forest...")
  log_msg("  {nrow(train_data)} échantillons × {length(feature_cols)} features", level = "info")

  # Préparer les données
  X_train <- train_data[, feature_cols, drop = FALSE]
  y_train <- as.factor(train_data[[target_col]])

  # Remplacement des NA par la médiane
  for (col in feature_cols) {
    na_mask <- is.na(X_train[[col]])
    if (any(na_mask)) {
      X_train[[col]][na_mask] <- median(X_train[[col]], na.rm = TRUE)
    }
  }

  # Retirer les colonnes à variance nulle
  var_check <- sapply(X_train, var, na.rm = TRUE)
  zero_var <- names(var_check[var_check == 0 | is.na(var_check)])
  if (length(zero_var) > 0) {
    log_msg("  Retrait de {length(zero_var)} features à variance nulle", level = "warning")
    X_train <- X_train[, !names(X_train) %in% zero_var, drop = FALSE]
    feature_cols <- setdiff(feature_cols, zero_var)
  }

  # Détermination automatique de mtry
  mtry <- CLASSIF_PARAMS$rf_mtry
  if (is.null(mtry)) {
    mtry <- floor(sqrt(ncol(X_train)))
  }

  log_msg("  ntree={CLASSIF_PARAMS$rf_ntree}, mtry={mtry}", level = "info")

  # Poids de classe pour compenser le déséquilibre
  case_wts <- compute_class_weights(y_train)

  # Entraînement avec ranger (plus rapide que randomForest)
  model <- ranger::ranger(
    x              = X_train,
    y              = y_train,
    num.trees      = CLASSIF_PARAMS$rf_ntree,
    mtry           = mtry,
    importance      = "impurity",
    probability     = TRUE,
    case.weights    = case_wts,
    seed            = CLASSIF_PARAMS$seed,
    verbose         = TRUE,
    num.threads     = parallel::detectCores() - 1
  )

  log_msg("Random Forest entraîné — OOB error : {round(model$prediction.error * 100, 2)}%",
          level = "success")

  # Stockage des métadonnées
  model$feature_cols <- feature_cols
  model$target_col   <- target_col
  model$n_classes    <- length(levels(y_train))
  model$class_names  <- levels(y_train)

  model
}

#' Prédiction avec le modèle Random Forest
#' @param model Modèle entraîné
#' @param test_data data.frame de test
#' @return Liste avec prédictions et probabilités
predict_rf <- function(model, test_data) {
  feature_cols <- model$feature_cols

  X_test <- test_data[, feature_cols, drop = FALSE]

  # Remplacement des NA
  for (col in feature_cols) {
    na_mask <- is.na(X_test[[col]])
    if (any(na_mask)) {
      X_test[[col]][na_mask] <- median(X_test[[col]], na.rm = TRUE)
    }
  }

  pred <- predict(model, data = X_test)

  list(
    predicted_class = pred$predictions |> apply(1, function(row) {
      model$class_names[which.max(row)]
    }),
    probabilities = pred$predictions
  )
}

# ==============================================================================
# 3. VALIDATION CROISÉE
# ==============================================================================

#' Validation croisée stratifiée k-fold
#' @param feature_matrix data.frame complet
#' @param feature_cols Colonnes features
#' @param target_col Colonne cible
#' @param k Nombre de folds
#' @param repeats Nombre de répétitions
#' @return Liste avec résultats par fold et résumé
cross_validate <- function(feature_matrix, feature_cols,
                            target_col = "species_name",
                            k = CLASSIF_PARAMS$cv_folds,
                            repeats = CLASSIF_PARAMS$cv_repeats) {
  log_msg("Validation croisée {k}-fold × {repeats} répétitions")

  set.seed(CLASSIF_PARAMS$seed)

  all_results <- list()

  for (rep in 1:repeats) {
    # Création des folds stratifiés
    species <- unique(feature_matrix[[target_col]])
    fold_assignments <- rep(NA, nrow(feature_matrix))

    for (sp in species) {
      sp_idx <- which(feature_matrix[[target_col]] == sp)
      n_sp <- length(sp_idx)
      fold_assignments[sp_idx] <- sample(rep(1:k, length.out = n_sp))
    }

    for (fold in 1:k) {
      log_msg("  Répétition {rep}/{repeats}, Fold {fold}/{k}")

      test_idx  <- which(fold_assignments == fold)
      train_idx <- which(fold_assignments != fold)

      train_fold <- feature_matrix[train_idx, ]
      test_fold  <- feature_matrix[test_idx, ]

      # Entraîner
      model <- train_random_forest(train_fold, feature_cols, target_col)

      # Prédire
      preds <- predict_rf(model, test_fold)

      # Évaluer
      y_true <- test_fold[[target_col]]
      y_pred <- preds$predicted_class

      accuracy <- mean(y_true == y_pred)

      all_results[[length(all_results) + 1]] <- data.frame(
        repeat_id = rep,
        fold      = fold,
        accuracy  = accuracy,
        n_test    = length(test_idx),
        stringsAsFactors = FALSE
      )
    }
  }

  results_df <- do.call(rbind, all_results)

  summary_stats <- list(
    mean_accuracy = mean(results_df$accuracy),
    sd_accuracy   = sd(results_df$accuracy),
    min_accuracy  = min(results_df$accuracy),
    max_accuracy  = max(results_df$accuracy),
    results       = results_df
  )

  log_msg(glue::glue(
    "CV terminée — Accuracy : {round(summary_stats$mean_accuracy*100, 2)}% ",
    "+/- {round(summary_stats$sd_accuracy*100, 2)}%"
  ), level = "success")

  summary_stats
}

# ==============================================================================
# 4. ÉVALUATION
# ==============================================================================

#' Calcul des métriques de classification détaillées
#' @param y_true Vecteur des vraies classes
#' @param y_pred Vecteur des classes prédites
#' @param class_names Noms des classes (optionnel)
#' @return Liste avec toutes les métriques
evaluate_classification <- function(y_true, y_pred, class_names = NULL) {
  log_msg("Évaluation de la classification...")

  # Matrice de confusion
  conf_mat <- table(Reference = y_true, Prediction = y_pred)

  if (is.null(class_names)) {
    class_names <- sort(unique(c(y_true, y_pred)))
  }

  # Overall Accuracy
  OA <- sum(diag(conf_mat)) / sum(conf_mat)

  # Métriques par classe
  n_classes <- length(class_names)
  per_class <- data.frame(
    species   = class_names,
    n_samples = as.integer(table(factor(y_true, levels = class_names))),
    stringsAsFactors = FALSE
  )

  per_class$precision <- sapply(class_names, function(cls) {
    if (!(cls %in% rownames(conf_mat) && cls %in% colnames(conf_mat))) return(0)
    tp <- conf_mat[cls, cls]
    fp <- sum(conf_mat[, cls]) - tp
    if (tp + fp == 0) 0 else tp / (tp + fp)
  })

  per_class$recall <- sapply(class_names, function(cls) {
    if (!(cls %in% rownames(conf_mat) && cls %in% colnames(conf_mat))) return(0)
    tp <- conf_mat[cls, cls]
    fn <- sum(conf_mat[cls, ]) - tp
    if (tp + fn == 0) 0 else tp / (tp + fn)
  })

  per_class$f1 <- with(per_class,
    ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  )

  # Kappa de Cohen — s'assurer que la matrice est carrée
  all_levels <- sort(unique(c(rownames(conf_mat), colnames(conf_mat))))
  full_conf <- matrix(0L, nrow = length(all_levels), ncol = length(all_levels),
                      dimnames = list(all_levels, all_levels))
  cr <- intersect(rownames(conf_mat), all_levels)
  cc <- intersect(colnames(conf_mat), all_levels)
  full_conf[cr, cc] <- conf_mat[cr, cc]

  p_o <- OA
  p_e <- sum(rowSums(full_conf) * colSums(full_conf)) / sum(full_conf)^2
  kappa <- (p_o - p_e) / (1 - p_e)

  # Macro-averaged metrics
  macro_precision <- mean(per_class$precision)
  macro_recall    <- mean(per_class$recall)
  macro_f1        <- mean(per_class$f1)

  # Weighted-averaged metrics
  weights <- per_class$n_samples / sum(per_class$n_samples)
  weighted_f1 <- sum(per_class$f1 * weights)

  results <- list(
    overall_accuracy  = OA,
    kappa             = kappa,
    macro_precision   = macro_precision,
    macro_recall      = macro_recall,
    macro_f1          = macro_f1,
    weighted_f1       = weighted_f1,
    per_class         = per_class,
    confusion_matrix  = conf_mat
  )

  # Affichage
  cli::cli_h2("Résultats de classification")
  cli::cli_alert_success("Overall Accuracy : {round(OA * 100, 2)}%")
  cli::cli_alert_info("Kappa de Cohen   : {round(kappa, 4)}")
  cli::cli_alert_info("Macro F1-Score   : {round(macro_f1 * 100, 2)}%")
  cli::cli_alert_info("Weighted F1      : {round(weighted_f1 * 100, 2)}%")

  results
}

#' Extraction de l'importance des variables depuis le modèle
#' @param model Modèle ranger
#' @return data.frame avec colonnes variable, importance
get_variable_importance <- function(model) {
  imp <- model$variable.importance

  data.frame(
    variable   = names(imp),
    importance = as.numeric(imp),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(importance))
}

# ==============================================================================
# 5. CNN TEMPOREL (optionnel, via torch)
# ==============================================================================

#' Définition d'un CNN 1D pour séries temporelles
#' Architecture : Conv1D → BatchNorm → ReLU → Conv1D → GAP → FC → Softmax
#' @param n_channels Nombre de canaux d'entrée (bandes spectrales)
#' @param n_timesteps Nombre de pas de temps
#' @param n_classes Nombre de classes de sortie
define_temporal_cnn <- function(n_channels = 10, n_timesteps = 73, n_classes = 21) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    cli::cli_alert_danger("Le package {.pkg torch} n'est pas installé.")
    cli::cli_text("Installer avec : {.code install.packages('torch'); torch::install_torch()}")
    return(NULL)
  }

  torch::nn_module(
    "TemporalCNN",

    initialize = function(n_channels, n_timesteps, n_classes) {
      # Bloc convolutif 1
      self$conv1 <- torch::nn_conv1d(n_channels, 64, kernel_size = 7, padding = 3)
      self$bn1   <- torch::nn_batch_norm1d(64)

      # Bloc convolutif 2
      self$conv2 <- torch::nn_conv1d(64, 128, kernel_size = 5, padding = 2)
      self$bn2   <- torch::nn_batch_norm1d(128)

      # Bloc convolutif 3
      self$conv3 <- torch::nn_conv1d(128, 256, kernel_size = 3, padding = 1)
      self$bn3   <- torch::nn_batch_norm1d(256)

      # Global Average Pooling + FC
      self$dropout <- torch::nn_dropout(p = 0.3)
      self$fc1 <- torch::nn_linear(256, 128)
      self$fc2 <- torch::nn_linear(128, n_classes)
    },

    forward = function(x) {
      # x shape: (batch, channels, timesteps)
      x <- self$conv1(x) |> self$bn1() |> torch::nnf_relu()
      x <- self$conv2(x) |> self$bn2() |> torch::nnf_relu()
      x <- self$conv3(x) |> self$bn3() |> torch::nnf_relu()

      # Global Average Pooling sur la dimension temporelle
      x <- x$mean(dim = 3)

      # Classification
      x <- self$dropout(x)
      x <- self$fc1(x) |> torch::nnf_relu()
      x <- self$dropout(x)
      x <- self$fc2(x)

      x
    }
  )
}

#' Préparation des données pour le CNN temporel
#' @param feature_matrix data.frame de features
#' @param feature_cols Colonnes features (séries temporelles des bandes)
#' @param target_col Colonne cible
#' @return Liste avec tenseurs X et Y
prepare_cnn_data <- function(feature_matrix, target_col = "species_name") {
  if (!requireNamespace("torch", quietly = TRUE)) {
    cli::cli_alert_danger("torch requis pour le CNN")
    return(NULL)
  }

  target_dates <- seq.Date(
    as.Date(TS_PARAMS$start_date),
    as.Date(TS_PARAMS$end_date),
    by = TS_PARAMS$target_interval_days
  )
  n_dates <- length(target_dates)
  doy_labels <- format(target_dates, "%j")

  # Extraire les séries temporelles par bande
  n_samples <- nrow(feature_matrix)
  n_bands   <- length(S2_BAND_NAMES)

  X_array <- array(0, dim = c(n_samples, n_bands, n_dates))

  for (b in seq_along(S2_BAND_NAMES)) {
    band <- S2_BAND_NAMES[b]
    col_names <- paste0(band, "_d", doy_labels)
    existing  <- col_names[col_names %in% names(feature_matrix)]

    if (length(existing) > 0) {
      mat <- as.matrix(feature_matrix[, existing])
      # Remplacement NA par 0
      mat[is.na(mat)] <- 0
      X_array[, b, 1:ncol(mat)] <- mat
    }
  }

  # Encodage des labels
  class_labels <- factor(feature_matrix[[target_col]])
  y_numeric <- as.integer(class_labels) - 1L  # 0-indexed

  list(
    X = torch::torch_tensor(X_array, dtype = torch::torch_float()),
    y = torch::torch_tensor(y_numeric, dtype = torch::torch_long()),
    class_names = levels(class_labels),
    n_channels  = n_bands,
    n_timesteps = n_dates,
    n_classes   = length(levels(class_labels))
  )
}

#' Entraînement du CNN temporel
#' @param train_tensors Sortie de prepare_cnn_data pour le train
#' @param val_tensors Sortie de prepare_cnn_data pour la validation
#' @param n_epochs Nombre d'époques
#' @param batch_size Taille de batch
#' @param lr Learning rate
#' @return Modèle entraîné
train_temporal_cnn <- function(train_tensors, val_tensors = NULL,
                                n_epochs = 50, batch_size = 32, lr = 0.001) {
  if (!requireNamespace("torch", quietly = TRUE)) return(NULL)

  log_msg("Entraînement CNN temporel...")
  log_msg("  Architecture : Conv1D(64) → Conv1D(128) → Conv1D(256) → GAP → FC",
          level = "info")
  log_msg("  Epochs={n_epochs}, Batch={batch_size}, LR={lr}", level = "info")

  # Créer le modèle
  model_fn <- define_temporal_cnn()
  model <- model_fn(
    n_channels  = train_tensors$n_channels,
    n_timesteps = train_tensors$n_timesteps,
    n_classes   = train_tensors$n_classes
  )

  optimizer <- torch::optim_adam(model$parameters, lr = lr)

  # Poids de classe pour compenser le déséquilibre dans la loss
  class_counts <- table(factor(
    as.integer(train_tensors$y$to(device = "cpu")) + 1L,
    levels = seq_len(train_tensors$n_classes)
  ))
  n_total <- sum(class_counts)
  n_cls   <- train_tensors$n_classes
  w_vec   <- n_total / (n_cls * pmax(class_counts, 1))
  class_weight_tensor <- torch::torch_tensor(as.numeric(w_vec), dtype = torch::torch_float())
  log_msg("  Poids de classe CNN : {paste(round(w_vec, 2), collapse=', ')}", level = "info")

  loss_fn <- torch::nn_cross_entropy_loss(weight = class_weight_tensor)

  # Dataset et DataLoader
  train_ds <- torch::tensor_dataset(train_tensors$X, train_tensors$y)
  train_dl <- torch::dataloader(train_ds, batch_size = batch_size, shuffle = TRUE)

  history <- data.frame(
    epoch      = integer(),
    train_loss = numeric(),
    train_acc  = numeric(),
    val_loss   = numeric(),
    val_acc    = numeric()
  )

  for (epoch in 1:n_epochs) {
    model$train()
    epoch_loss <- 0
    epoch_correct <- 0
    epoch_total <- 0

    coro::loop(for (batch in train_dl) {
      optimizer$zero_grad()

      output <- model(batch[[1]])
      loss <- loss_fn(output, batch[[2]])

      loss$backward()
      optimizer$step()

      epoch_loss <- epoch_loss + loss$item()
      preds <- output$argmax(dim = 2)
      epoch_correct <- epoch_correct + (preds == batch[[2]])$sum()$item()
      epoch_total <- epoch_total + length(batch[[2]])
    })

    train_acc <- epoch_correct / epoch_total
    avg_loss  <- epoch_loss / length(train_dl)

    # Validation
    val_loss <- NA
    val_acc  <- NA
    if (!is.null(val_tensors)) {
      model$eval()
      torch::with_no_grad({
        val_output <- model(val_tensors$X)
        val_loss <- loss_fn(val_output, val_tensors$y)$item()
        val_preds <- val_output$argmax(dim = 2)
        val_acc <- (val_preds == val_tensors$y)$sum()$item() / length(val_tensors$y)
      })
    }

    history <- rbind(history, data.frame(
      epoch = epoch, train_loss = avg_loss, train_acc = train_acc,
      val_loss = val_loss, val_acc = val_acc
    ))

    if (epoch %% 10 == 0 || epoch == 1) {
      msg <- glue::glue("  Epoch {epoch}/{n_epochs} — Loss: {round(avg_loss, 4)}, Acc: {round(train_acc*100, 1)}%")
      if (!is.na(val_acc)) {
        msg <- paste0(msg, glue::glue(" | Val Acc: {round(val_acc*100, 1)}%"))
      }
      log_msg(msg)
    }
  }

  log_msg("CNN entraîné — Accuracy finale : {round(tail(history$train_acc,1)*100, 1)}%",
          level = "success")

  list(model = model, history = history, class_names = train_tensors$class_names)
}

# ==============================================================================
# 6. SAUVEGARDE / CHARGEMENT DES MODÈLES
# ==============================================================================

#' Sauvegarde du modèle et des métadonnées
#' @param model Modèle (ranger ou torch)
#' @param eval_results Résultats d'évaluation
#' @param model_name Nom du modèle
#' @param output_dir Répertoire de sortie
save_model <- function(model, eval_results = NULL, model_name = "treesatai_rf",
                        output_dir = NULL) {
  if (is.null(output_dir)) output_dir <- file.path(.get_project_root(), "output", "models")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Sauvegarder le modèle R
  model_path <- file.path(output_dir, paste0(model_name, ".rds"))
  saveRDS(model, model_path)

  # Sauvegarder les résultats d'évaluation
  if (!is.null(eval_results)) {
    eval_path <- file.path(output_dir, paste0(model_name, "_evaluation.rds"))
    saveRDS(eval_results, eval_path)

    # CSV des métriques par classe
    csv_path <- file.path(output_dir, paste0(model_name, "_per_class_metrics.csv"))
    readr::write_csv(eval_results$per_class, csv_path)
  }

  # Métadonnées JSON
  metadata <- list(
    model_name    = model_name,
    date_created  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_classes     = if (!is.null(model$n_classes)) model$n_classes else NA,
    n_features    = if (!is.null(model$feature_cols)) length(model$feature_cols) else NA,
    oob_error     = if (!is.null(model$prediction.error)) model$prediction.error else NA,
    overall_accuracy = if (!is.null(eval_results)) eval_results$overall_accuracy else NA,
    kappa         = if (!is.null(eval_results)) eval_results$kappa else NA,
    class_names   = if (!is.null(model$class_names)) model$class_names else SPECIES$latin,
    class_weighting = TRUE
  )

  meta_path <- file.path(output_dir, paste0(model_name, "_metadata.json"))
  jsonlite::write_json(metadata, meta_path, pretty = TRUE, auto_unbox = TRUE)

  log_msg("Modèle sauvegardé : {model_path}", level = "success")
  invisible(model_path)
}

#' Chargement d'un modèle sauvegardé
#' @param model_path Chemin vers le fichier .rds
#' @return Modèle chargé
load_model <- function(model_path) {
  if (!file.exists(model_path)) {
    cli::cli_alert_danger("Fichier modèle introuvable : {model_path}")
    return(NULL)
  }
  model <- readRDS(model_path)
  log_msg("Modèle chargé : {model_path}", level = "success")
  model
}

