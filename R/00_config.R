# ==============================================================================
# TreeSatAI-Time-Series (IGNF) — Configuration
# Classification d'essences forestières par séries temporelles Sentinel-2
# 20 espèces européennes identifiées par signature phénologique
# ==============================================================================

# --- Chemins du projet -------------------------------------------------------
# En mode package, les chemins sont relatifs au working directory
# L'utilisateur peut les changer via options(treesatnemeton.data_dir = "...")
.get_project_root <- function() {
  getOption("treesatnemeton.project_root",
            default = if (requireNamespace("here", quietly = TRUE)) here::here() else getwd())
}

PROJECT_ROOT  <- .get_project_root()
DATA_DIR      <- file.path(PROJECT_ROOT, "data")
RAW_DIR       <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
TS_DIR        <- file.path(DATA_DIR, "timeseries")
OUTPUT_DIR    <- file.path(PROJECT_ROOT, "output")
FIGURES_DIR   <- file.path(PROJECT_ROOT, "figures")
MODELS_DIR    <- file.path(OUTPUT_DIR, "models")

#' Initialise les répertoires du projet
#' @param root Répertoire racine (par défaut : working directory)
#' @export
init_project_dirs <- function(root = .get_project_root()) {
  dirs <- c(
    file.path(root, "data"),
    file.path(root, "data", "raw"),
    file.path(root, "data", "processed"),
    file.path(root, "data", "timeseries"),
    file.path(root, "output"),
    file.path(root, "output", "models"),
    file.path(root, "figures")
  )
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  invisible(dirs)
}

# --- Paramètres Sentinel-2 ---------------------------------------------------
S2_BANDS <- list(
  B02 = list(name = "Blue",     wavelength = 490,  resolution = 10),
  B03 = list(name = "Green",    wavelength = 560,  resolution = 10),
  B04 = list(name = "Red",      wavelength = 665,  resolution = 10),
  B05 = list(name = "RedEdge1", wavelength = 705,  resolution = 20),
  B06 = list(name = "RedEdge2", wavelength = 740,  resolution = 20),
  B07 = list(name = "RedEdge3", wavelength = 783,  resolution = 20),
  B08 = list(name = "NIR",      wavelength = 842,  resolution = 10),
  B8A = list(name = "NIR2",     wavelength = 865,  resolution = 20),
  B11 = list(name = "SWIR1",    wavelength = 1610, resolution = 20),
  B12 = list(name = "SWIR2",    wavelength = 2190, resolution = 20)
)

# Bandes utilisées pour la classification (10 bandes spectrales)
S2_BAND_NAMES <- names(S2_BANDS)

# --- Les 20 espèces européennes TreeSatAI ------------------------------------
SPECIES <- data.frame(
  code = 1:20,
  latin = c(
    "Quercus robur",           # Chêne pédonculé
    "Quercus petraea",         # Chêne sessile
    "Quercus pubescens",       # Chêne pubescent
    "Quercus ilex",            # Chêne vert
    "Fagus sylvatica",         # Hêtre
    "Castanea sativa",         # Châtaignier
    "Carpinus betulus",        # Charme
    "Betula pendula",          # Bouleau verruqueux
    "Fraxinus excelsior",      # Frêne commun
    "Acer pseudoplatanus",     # Érable sycomore
    "Populus spp.",             # Peupliers
    "Robinia pseudoacacia",    # Robinier faux-acacia
    "Picea abies",             # Épicéa commun
    "Abies alba",              # Sapin pectiné
    "Pseudotsuga menziesii",   # Douglas
    "Pinus sylvestris",        # Pin sylvestre
    "Pinus pinaster",          # Pin maritime
    "Pinus nigra",             # Pin noir
    "Pinus halepensis",        # Pin d'Alep
    "Larix decidua"            # Mélèze d'Europe
  ),
  french = c(
    "Chêne pédonculé", "Chêne sessile", "Chêne pubescent",
    "Chêne vert", "Hêtre", "Châtaignier", "Charme",
    "Bouleau verruqueux", "Frêne commun", "Érable sycomore",
    "Peupliers", "Robinier faux-acacia", "Épicéa commun",
    "Sapin pectiné", "Douglas", "Pin sylvestre",
    "Pin maritime", "Pin noir", "Pin d'Alep", "Mélèze d'Europe"
  ),
  type = c(
    rep("feuillu", 12),
    rep("résineux", 8)
  ),
  phenologie = c(
    rep("caducifolié", 3),  # Chênes caducs
    "sempervirent",          # Chêne vert
    rep("caducifolié", 8),  # Autres feuillus caducs
    rep("sempervirent", 6), # Résineux persistants
    "sempervirent",          # Pin d'Alep
    "caducifolié"            # Mélèze (résineux caduc !)
  ),
  stringsAsFactors = FALSE
)

# Palette de couleurs par espèce
SPECIES_COLORS <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3",
  "#bebada", "#fb8072", "#003c30", "#01665e", "#35978f",
  "#80cdc1", "#c7eae5", "#543005", "#bf812d", "#dfc27d"
)
names(SPECIES_COLORS) <- SPECIES$latin

# --- Paramètres de la série temporelle ----------------------------------------
TS_PARAMS <- list(
  # Période d'acquisition (1 année complète)
  year          = 2021,
  start_date    = "2021-01-01",
  end_date      = "2021-12-31",

  # Résolution temporelle cible après interpolation
  target_interval_days = 5,   # Interpolation à 5 jours
  n_dates_target       = 73,  # ~365/5 = 73 dates

  # Seuil de couverture nuageuse max par scène
  max_cloud_cover = 30,       # en %

  # Filtre Savitzky-Golay pour lissage
  sg_filter_order  = 3,       # Ordre du polynôme
  sg_filter_length = 7        # Taille de la fenêtre (impair)
)

# --- Indices spectraux à calculer --------------------------------------------
SPECTRAL_INDICES <- c("NDVI", "EVI", "NDWI", "CRI", "RENDVI", "NBR")

# --- Paramètres de classification --------------------------------------------
CLASSIF_PARAMS <- list(
  # Split train/test
  train_ratio = 0.7,
  seed        = 42,

  # Random Forest
  rf_ntree     = 500,
  rf_mtry      = NULL,  # auto : sqrt(n_features)
  rf_importance = TRUE,

  # Validation croisée
  cv_folds     = 5,
  cv_repeats   = 3,

  # Features à utiliser
  use_spectral_bands  = TRUE,  # Bandes brutes interpolées
 use_spectral_indices = TRUE,  # Indices spectraux calculés
  use_phenometrics     = TRUE,  # Métriques phénologiques extraites
  use_temporal_stats   = TRUE   # Statistiques temporelles (mean, sd, etc.)
)

# --- Paramètres de visualisation ---------------------------------------------
VIS_PARAMS <- list(
  dpi        = 300,
  width_cm   = 20,
  height_cm  = 15,
  theme_base = "theme_minimal",
  font_size  = 10
)

# --- Message de bienvenue (affiché via .onAttach dans zzz.R) ----------------
