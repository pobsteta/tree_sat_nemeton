#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Pont R ↔ Python (reticulate)
#
# Interface R pour utiliser les modèles PyTorch via reticulate.
# Pattern identique à flairhub et open_canopy.
#
# Usage :
#   source("R/09_python_bridge.R")
#   setup_python_env()   # Vérifie/installe le conda env "treesat"
#   model <- py_load_model("output/models/treesatai_tempcnn_best.pt")
#   preds <- py_predict_pixels(model, pixel_array)
# ==============================================================================

# Configuration chargée via le package

# ==============================================================================
# 1. CONFIGURATION DE L'ENVIRONNEMENT PYTHON
# ==============================================================================

#' Initialiser l'environnement Python conda "treesat"
#'
#' Vérifie que l'env conda existe, le crée si besoin à partir de environment.yml,
#' puis configure reticulate pour l'utiliser.
#'
#' @param conda_env Nom de l'environnement conda (défaut: "treesat")
#' @param force_reinstall Re-créer l'environnement même s'il existe
#' @return TRUE si OK, FALSE sinon
setup_python_env <- function(conda_env = "treesat", force_reinstall = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Le package 'reticulate' est requis. Installez-le : install.packages('reticulate')")
  }

  cli::cli_h2("Configuration Python ({conda_env})")

  # Vérifier si conda est disponible
  conda_bin <- tryCatch(
    reticulate::conda_binary(),
    error = function(e) NULL
  )

  if (is.null(conda_bin)) {
    cli::cli_alert_danger("Conda non trouvé.")
    cli::cli_text("Installez Miniconda : {.code reticulate::install_miniconda()}")
    cli::cli_text("Ou installez Anaconda/Miniconda manuellement.")
    return(FALSE)
  }

  cli::cli_alert_info("Conda : {conda_bin}")

  # Lister les envs existants
  envs <- reticulate::conda_list()
  env_exists <- conda_env %in% envs$name

  if (env_exists && !force_reinstall) {
    cli::cli_alert_success("Environnement '{conda_env}' trouvé")
  } else {
    # Créer l'environnement depuis environment.yml
    yml_path <- file.path(here::here(), "environment.yml")

    if (file.exists(yml_path)) {
      cli::cli_alert_info("Création de l'environnement depuis environment.yml...")
      cli::cli_text("(cela peut prendre plusieurs minutes)")

      tryCatch({
        system2(
          conda_bin,
          args = c("env", "create", "-f", shQuote(yml_path),
                   if (force_reinstall) "--force" else NULL),
          stdout = TRUE, stderr = TRUE
        )
        cli::cli_alert_success("Environnement '{conda_env}' créé")
      }, error = function(e) {
        cli::cli_alert_danger("Erreur lors de la création de l'env : {e$message}")
        cli::cli_text("Essayez manuellement : conda env create -f environment.yml")
        return(FALSE)
      })
    } else {
      # Créer un env minimal sans le yml
      cli::cli_alert_warning("environment.yml non trouvé — création d'un env minimal")
      reticulate::conda_create(
        envname = conda_env,
        packages = c("python=3.11", "numpy", "pandas", "scikit-learn")
      )
      # Installer PyTorch et dépendances via pip
      reticulate::conda_install(
        envname = conda_env,
        packages = c("torch", "torchvision", "onnx", "onnxruntime"),
        pip = TRUE
      )
    }
  }

  # Activer l'environnement
  reticulate::use_condaenv(conda_env, required = TRUE)
  cli::cli_alert_success("Python activé : {reticulate::py_config()$python}")

  # Vérifier les imports critiques
  imports_ok <- tryCatch({
    reticulate::py_module_available("torch") &&
      reticulate::py_module_available("numpy")
  }, error = function(e) FALSE)

  if (!imports_ok) {
    cli::cli_alert_danger("PyTorch ou NumPy non disponible dans l'env '{conda_env}'")
    cli::cli_text("Réinstallez : {.code setup_python_env(force_reinstall = TRUE)}")
    return(FALSE)
  }

  torch_version <- reticulate::py_eval("__import__('torch').__version__")
  cuda_available <- reticulate::py_eval("__import__('torch').cuda.is_available()")
  device <- if (cuda_available) "GPU (CUDA)" else "CPU"

  cli::cli_alert_success("PyTorch {torch_version} — {device}")

  invisible(TRUE)
}


# ==============================================================================
# 2. CHARGEMENT ET PRÉDICTION AVEC LES MODÈLES PYTORCH
# ==============================================================================

#' Charger un modèle PyTorch depuis un checkpoint .pt
#'
#' @param model_path Chemin vers le fichier .pt (sauvegardé par train.py)
#' @param device "cuda" ou "cpu" (NULL = auto)
#' @return Objet Python (dict) contenant le modèle, les noms de classes, etc.
py_load_model <- function(model_path, device = NULL) {
  if (!file.exists(model_path)) {
    stop("Modèle non trouvé : ", model_path)
  }

  # Importer le module predict de Python
  py_predict <- import_predict_module()

  if (is.null(device)) {
    device <- reticulate::py_eval(
      "str(__import__('torch').device('cuda' if __import__('torch').cuda.is_available() else 'cpu'))"
    )
  }

  cli::cli_alert_info("Chargement du modèle : {basename(model_path)}")

  result <- py_predict$load_model(model_path, device)
  model <- result[[1]]
  class_names <- result[[2]]
  checkpoint <- result[[3]]

  model_name <- tryCatch(
    as.character(checkpoint$get("model_name", "unknown")),
    error = function(e) "unknown"
  )
  val_acc <- tryCatch(
    as.numeric(checkpoint$get("val_acc", 0)),
    error = function(e) NA_real_
  )

  cli::cli_alert_success("Modèle '{model_name}' chargé (val_acc={round(val_acc * 100, 1)}%)")

  list(
    model       = model,
    class_names = as.character(reticulate::py_to_r(class_names)),
    checkpoint  = checkpoint,
    model_name  = model_name,
    device      = device
  )
}


#' Prédiction PyTorch pixel par pixel
#'
#' Convertit les données R en tensors NumPy, appelle predict_pixels() de Python,
#' et renvoie les résultats en R.
#'
#' @param py_model Objet retourné par py_load_model()
#' @param pixel_data Matrice R (n_pixels × n_bands × n_timesteps) ou array 3D
#' @param batch_size Taille de batch pour l'inférence
#' @return Liste R avec predicted_class, probabilities, max_proba
py_predict_pixels <- function(py_model, pixel_data, batch_size = 512L) {
  py_predict <- import_predict_module()
  np <- reticulate::import("numpy", convert = FALSE)

  # Convertir en numpy array
  if (is.matrix(pixel_data)) {
    # Si c'est une matrice 2D (n_pixels × features), on la reshape
    cli::cli_alert_warning("Données 2D détectées — le modèle attend (n_pixels, n_bands, n_timesteps)")
    np_data <- np$array(pixel_data, dtype = "float32")
  } else if (is.array(pixel_data) && length(dim(pixel_data)) == 3) {
    np_data <- np$array(pixel_data, dtype = "float32")
  } else {
    np_data <- np$array(as.array(pixel_data), dtype = "float32")
  }

  cli::cli_alert_info("Prédiction PyTorch sur {dim(pixel_data)[1]} pixels...")

  results <- py_predict$predict_pixels(
    model      = py_model$model,
    pixel_data = np_data,
    batch_size = as.integer(batch_size),
    device     = NULL
  )

  # Convertir en R
  list(
    predicted_class = as.integer(reticulate::py_to_r(results$predicted_class)) + 1L,  # 0-based → 1-based
    probabilities   = reticulate::py_to_r(results$probabilities),
    max_proba       = as.numeric(reticulate::py_to_r(results$max_proba))
  )
}


#' Prédiction multi-source (S2 + S1) avec PyTorch
#'
#' @param py_model Objet retourné par py_load_model()
#' @param s2_data Array 3D (n_pixels, n_s2_bands, n_s2_timesteps)
#' @param s1_data Array 3D (n_pixels, n_s1_bands, n_s1_timesteps) ou NULL
#' @param batch_size Taille de batch
#' @return Liste R avec predicted_class, probabilities, max_proba
py_predict_multisource <- function(py_model, s2_data, s1_data = NULL, batch_size = 512L) {
  py_predict <- import_predict_module()
  np <- reticulate::import("numpy", convert = FALSE)

  np_s2 <- np$array(s2_data, dtype = "float32")
  np_s1 <- if (!is.null(s1_data)) np$array(s1_data, dtype = "float32") else NULL

  cli::cli_alert_info("Prédiction multi-source PyTorch sur {dim(s2_data)[1]} pixels...")

  results <- py_predict$predict_multisource(
    model      = py_model$model,
    s2_data    = np_s2,
    s1_data    = np_s1,
    batch_size = as.integer(batch_size),
    device     = NULL
  )

  list(
    predicted_class = as.integer(reticulate::py_to_r(results$predicted_class)) + 1L,
    probabilities   = reticulate::py_to_r(results$probabilities),
    max_proba       = as.numeric(reticulate::py_to_r(results$max_proba))
  )
}


# ==============================================================================
# 3. ENTRAÎNEMENT DEPUIS R
# ==============================================================================

#' Lancer l'entraînement PyTorch depuis R
#'
#' Appelle python/train.py avec les arguments spécifiés.
#'
#' @param data_path Chemin vers feature_matrix.csv
#' @param model_type Architecture : "tempcnn", "lstm", "transformer", "inception", "multisource"
#' @param epochs Nombre d'époques
#' @param batch_size Taille de batch
#' @param lr Learning rate
#' @param patience Early stopping patience
#' @param conda_env Nom de l'env conda
#' @return Chemin vers le meilleur modèle sauvegardé
py_train_model <- function(data_path,
                           model_type = "tempcnn",
                           epochs     = 100L,
                           batch_size = 64L,
                           lr         = 1e-3,
                           patience   = 15L,
                           conda_env  = "treesat") {
  if (!file.exists(data_path)) {
    stop("Données non trouvées : ", data_path)
  }

  cli::cli_h2("Entraînement PyTorch — {model_type}")

  train_script <- file.path(here::here(), "python", "train.py")
  if (!file.exists(train_script)) {
    stop("Script d'entraînement non trouvé : ", train_script)
  }

  # Trouver le python du conda env
  envs <- reticulate::conda_list()
  env_info <- envs[envs$name == conda_env, ]
  if (nrow(env_info) == 0) {
    stop("Environnement conda '", conda_env, "' non trouvé. Lancez setup_python_env() d'abord.")
  }
  python_bin <- env_info$python[1]

  # Construire la commande
  cmd <- sprintf(
    '%s "%s" --data "%s" --model %s --epochs %d --batch-size %d --lr %s --patience %d',
    shQuote(python_bin), train_script, data_path,
    model_type, as.integer(epochs), as.integer(batch_size),
    format(lr, scientific = FALSE), as.integer(patience)
  )

  cli::cli_alert_info("Commande : {cmd}")
  cli::cli_text("")

  # Exécuter
  t_start <- Sys.time()
  exit_code <- system(cmd)
  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

  if (exit_code != 0) {
    cli::cli_alert_danger("Entraînement échoué (code {exit_code})")
    return(NULL)
  }

  # Trouver le modèle sauvegardé
  model_path <- file.path(here::here(), "output", "models",
                           paste0("treesatai_", model_type, "_best.pt"))

  if (file.exists(model_path)) {
    cli::cli_alert_success("Modèle sauvegardé : {model_path}")
    cli::cli_alert_info("Durée : {round(t_elapsed, 1)} minutes")
  } else {
    cli::cli_alert_warning("Modèle attendu non trouvé : {model_path}")
  }

  model_path
}


# ==============================================================================
# 4. EXPORT ONNX DEPUIS R
# ==============================================================================

#' Exporter un modèle PyTorch en ONNX
#'
#' @param model_path Chemin vers le .pt
#' @param output_path Chemin de sortie .onnx (NULL = même nom que le .pt)
#' @param n_channels Nombre de bandes spectrales
#' @param n_timesteps Nombre de timesteps
#' @return Chemin du fichier ONNX
py_export_onnx <- function(model_path, output_path = NULL,
                           n_channels = 10L, n_timesteps = 73L) {
  py_predict <- import_predict_module()

  if (is.null(output_path)) {
    output_path <- sub("\\.pt$", ".onnx", model_path)
  }

  cli::cli_alert_info("Export ONNX : {basename(model_path)} → {basename(output_path)}")

  py_model <- py_load_model(model_path)
  torch <- reticulate::import("torch")

  model_cpu <- py_model$model$cpu()
  py_predict$export_onnx(
    model_cpu, output_path,
    n_channels  = as.integer(n_channels),
    n_timesteps = as.integer(n_timesteps)
  )

  cli::cli_alert_success("ONNX exporté : {output_path}")
  output_path
}


# ==============================================================================
# 5. CONVERSION DES FEATURES R → FORMAT PYTORCH
# ==============================================================================

#' Convertir une feature matrix R (pixels × features) en format PyTorch
#'
#' Le modèle Python attend (n_pixels, n_bands, n_timesteps).
#' Cette fonction reconstruit l'array 3D à partir de la matrice plate.
#'
#' @param feature_matrix Matrice de features (sortie de extract_pixel_features)
#' @param bands Noms des bandes à utiliser
#' @param n_timesteps Nombre de timesteps par bande
#' @return Array 3D (n_pixels, n_bands, n_timesteps) prêt pour PyTorch
features_to_pytorch_format <- function(feature_matrix, bands = NULL, n_timesteps = 73L) {
  if (is.null(bands)) {
    bands <- S2_BAND_NAMES
  }

  n_pixels <- nrow(feature_matrix)
  n_bands  <- length(bands)

  cli::cli_alert_info("Conversion R → PyTorch : {n_pixels} pixels × {n_bands} bandes × {n_timesteps} timesteps")

  # Chercher les colonnes de séries temporelles brutes
  # On s'attend à des colonnes comme "B02_t001", "B02_t002", ... "B02_t073"
  ts_array <- array(0, dim = c(n_pixels, n_bands, n_timesteps))

  for (b in seq_along(bands)) {
    band <- bands[b]
    # Chercher les colonnes de cette bande (format ts_xxx ou t_xxx)
    col_patterns <- c(
      paste0(band, "_t", sprintf("%03d", 1:n_timesteps)),  # B02_t001
      paste0(band, "_ts_", 1:n_timesteps),                 # B02_ts_1
      paste0(band, "_", 1:n_timesteps)                     # B02_1
    )

    # Essayer chaque pattern
    found <- FALSE
    for (pattern_set in list(
      paste0(band, "_t", sprintf("%03d", 1:n_timesteps)),
      paste0(band, "_ts_", 1:n_timesteps),
      paste0(band, "_", 1:n_timesteps)
    )) {
      matched <- pattern_set %in% colnames(feature_matrix)
      if (sum(matched) >= n_timesteps * 0.5) {
        for (t in seq_len(n_timesteps)) {
          if (pattern_set[t] %in% colnames(feature_matrix)) {
            ts_array[, b, t] <- feature_matrix[, pattern_set[t]]
          }
        }
        found <- TRUE
        break
      }
    }

    if (!found) {
      cli::cli_alert_warning("Bande '{band}' : séries temporelles non trouvées dans les features")
    }
  }

  # Remplacer NA par 0
  ts_array[is.na(ts_array)] <- 0

  ts_array
}


#' Convertir un cube S2 (sortie de build_s2_cube) en array PyTorch
#'
#' @param cube_list Liste de SpatRaster (1 par date)
#' @param valid_idx Indices des pixels valides
#' @param bands Bandes à extraire
#' @return Array 3D (n_pixels, n_bands, n_timesteps)
cube_to_pytorch_array <- function(cube_list, valid_idx, bands = NULL) {
  if (is.null(bands)) {
    bands <- S2_BAND_NAMES
  }

  n_dates  <- length(cube_list)
  n_pixels <- length(valid_idx)
  n_bands  <- length(bands)

  cli::cli_alert_info("Extraction cube → array : {n_pixels} pixels × {n_bands} bandes × {n_dates} dates")

  ts_array <- array(0, dim = c(n_pixels, n_bands, n_dates))

  for (b in seq_along(bands)) {
    band <- bands[b]
    for (d in seq_len(n_dates)) {
      layer_name <- grep(paste0("^", band, "_"), names(cube_list[[d]]), value = TRUE)
      if (length(layer_name) > 0) {
        all_vals <- terra::values(cube_list[[d]][[layer_name[1]]])
        ts_array[, b, d] <- as.numeric(all_vals[valid_idx])
      }
    }
  }

  # Remplacer NA par 0
  ts_array[is.na(ts_array)] <- 0

  ts_array
}


# ==============================================================================
# 6. CLASSIFICATION PYTORCH INTÉGRÉE AU PIPELINE R
# ==============================================================================

#' Classification des pixels via PyTorch (alternative à classify_pixels avec ranger)
#'
#' @param pixel_features Sortie de extract_pixel_features() ou cube + valid_idx
#' @param model_path Chemin vers le modèle .pt
#' @param cube_list Cube S2 (si pixel_features ne contient pas les TS brutes)
#' @param batch_size Taille de batch pour l'inférence
#' @return Liste compatible avec le format de classify_pixels()
classify_pixels_pytorch <- function(pixel_features, model_path,
                                     cube_list = NULL, batch_size = 512L) {
  cli::cli_h3("Classification PyTorch")

  # Charger le modèle
  py_model <- py_load_model(model_path)

  # Préparer les données
  if (!is.null(cube_list)) {
    # Extraire directement du cube (séries temporelles brutes)
    pixel_data <- cube_to_pytorch_array(
      cube_list, pixel_features$valid_idx, bands = S2_BAND_NAMES
    )
  } else {
    # Convertir depuis la matrice de features
    pixel_data <- features_to_pytorch_format(pixel_features$features)
  }

  cli::cli_alert_info("Shape des données : {paste(dim(pixel_data), collapse = ' × ')}")

  # Prédiction
  results <- py_predict_pixels(py_model, pixel_data, batch_size = batch_size)

  # Formater en sortie compatible avec classify_pixels()
  list(
    class_idx    = results$predicted_class,
    class_names  = py_model$class_names,
    max_proba    = results$max_proba,
    all_probas   = results$probabilities,
    valid_idx    = pixel_features$valid_idx
  )
}


# ==============================================================================
# 7. UTILITAIRES INTERNES
# ==============================================================================

#' Importer le module python/predict.py via reticulate
#' @return Module Python
import_predict_module <- function() {
  python_dir <- file.path(here::here(), "python")

  # Ajouter le répertoire python au sys.path
  sys <- reticulate::import("sys")
  if (!(python_dir %in% sys$path)) {
    sys$path$insert(0L, python_dir)
  }

  # Importer le module predict
  reticulate::import_from_path("predict", path = python_dir)
}


#' Vérifier que l'environnement Python est prêt
#' @return TRUE si OK
check_python_ready <- function() {
  tryCatch({
    reticulate::py_module_available("torch")
  }, error = function(e) {
    cli::cli_alert_danger("Python non configuré. Lancez d'abord : setup_python_env()")
    FALSE
  })
}


