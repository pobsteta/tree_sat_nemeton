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

# Facteur d'échelle Sentinel-2 L2A (réflectance entière → 0-1)
S2_SCALE_FACTOR <- 10000

# Bandes utilisées pour la classification (10 bandes spectrales)
S2_BAND_NAMES <- names(S2_BANDS)

# --- Les 20 espèces européennes TreeSatAI ------------------------------------
SPECIES <- data.frame(
  code = 1:21,
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
    "Larix decidua",           # Mélèze d'Europe
    "Cleared"                  # Coupe rase / vide forestier
  ),
  french = c(
    "Chêne pédonculé", "Chêne sessile", "Chêne pubescent",
    "Chêne vert", "Hêtre", "Châtaignier", "Charme",
    "Bouleau verruqueux", "Frêne commun", "Érable sycomore",
    "Peupliers", "Robinier faux-acacia", "Épicéa commun",
    "Sapin pectiné", "Douglas", "Pin sylvestre",
    "Pin maritime", "Pin noir", "Pin d'Alep", "Mélèze d'Europe",
    "Coupe/Vide"
  ),
  type = c(
    rep("feuillu", 12),
    rep("résineux", 8),
    "non-boisé"                # Cleared : ni feuillu ni résineux
  ),
  phenologie = c(
    rep("caducifolié", 3),  # Chênes caducs
    "sempervirent",          # Chêne vert
    rep("caducifolié", 8),  # Autres feuillus caducs
    rep("sempervirent", 6), # Résineux persistants
    "sempervirent",          # Pin d'Alep
    "caducifolié",           # Mélèze (résineux caduc !)
    NA_character_            # Cleared : pas de phénologie
  ),
  stringsAsFactors = FALSE
)

# --- Regroupement en 10 classes forestières ------------------------------------
# Mapping : nom français (20 espèces) → nom du groupe (10 classes)
SPECIES_GROUPS <- c(
  "Chêne pédonculé"      = "Chênes caducs",
  "Chêne sessile"        = "Chênes caducs",
  "Chêne pubescent"      = "Chênes caducs",
  "Chêne vert"           = "Chêne vert",
  "Hêtre"                = "Hêtre",
  "Châtaignier"          = "Châtaignier",
  "Charme"               = "Autres feuillus",
  "Bouleau verruqueux"   = "Autres feuillus",
  "Frêne commun"         = "Autres feuillus",
  "Érable sycomore"      = "Autres feuillus",
  "Peupliers"            = "Autres feuillus",
  "Robinier faux-acacia" = "Autres feuillus",
  "Épicéa commun"        = "Sapins-Épicéas",
  "Sapin pectiné"        = "Sapins-Épicéas",
  "Douglas"              = "Douglas",
  "Pin sylvestre"        = "Pin sylvestre",
  "Pin maritime"         = "Autres Pins",
  "Pin noir"             = "Autres Pins",
  "Pin d'Alep"           = "Autres Pins",
  "Mélèze d'Europe"      = "Mélèze",
  "Coupe/Vide"           = "Coupe/Vide"
)

# Metadata des 10 groupes
SPECIES_GROUPS_INFO <- data.frame(
  group = c("Chênes caducs", "Chêne vert", "Hêtre", "Châtaignier",
            "Autres feuillus", "Sapins-Épicéas", "Douglas",
            "Pin sylvestre", "Autres Pins", "Mélèze", "Coupe/Vide"),
  type = c("feuillu", "feuillu", "feuillu", "feuillu", "feuillu",
           "résineux", "résineux", "résineux", "résineux", "résineux",
           "non-boisé"),
  phenologie = c("caducifolié", "sempervirent", "caducifolié", "caducifolié",
                 "caducifolié", "sempervirent", "sempervirent",
                 "sempervirent", "sempervirent", "caducifolié",
                 NA_character_),
  stringsAsFactors = FALSE
)

# Palette de couleurs par groupe (10 classes)
SPECIES_GROUP_COLORS <- c(
  "Chênes caducs"   = "#1b9e77",
  "Chêne vert"      = "#d95f02",
  "Hêtre"           = "#7570b3",
  "Châtaignier"     = "#e7298a",
  "Autres feuillus" = "#66a61e",
  "Sapins-Épicéas"  = "#003c30",
  "Douglas"         = "#01665e",
  "Pin sylvestre"   = "#80cdc1",
  "Autres Pins"     = "#bf812d",
  "Mélèze"          = "#dfc27d",
  "Coupe/Vide"      = "#969696"
)

# Palette de couleurs par espèce (21 classes : 20 essences + Cleared)
SPECIES_COLORS <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3",
  "#bebada", "#fb8072", "#003c30", "#01665e", "#35978f",
  "#80cdc1", "#c7eae5", "#543005", "#bf812d", "#dfc27d",
  "#969696"   # Cleared — gris
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
  use_spectral_bands   = TRUE,  # Bandes brutes interpolées
  use_spectral_indices = TRUE,  # Indices spectraux calculés
  use_phenometrics     = TRUE,  # Métriques phénologiques extraites
  use_temporal_stats   = TRUE,  # Statistiques temporelles (mean, sd, etc.)
  use_terrain          = TRUE,  # Features topographiques (MNT, pente, exposition, TWI)

  # --- Classification soft (probabilités) ---
  # Exporter les probabilités complètes par essence (raster multi-bandes)
  export_probabilities = TRUE,

  # Seuil de présence : une essence est "présente" si sa proba dépasse ce seuil
  presence_threshold = 0.10,

  # Exporter la carte d'entropie de Shannon (taux de mélange)
  export_shannon = TRUE,

  # Exporter les cartes de présence binaires par essence (proba > seuil)
  export_presence = TRUE
)

# --- Paramètres MNT / terrain -------------------------------------------------
DEM_PARAMS <- list(
  # Source DEM pour la prédiction spatiale
  # "copernicus" = Copernicus DEM 30 m (Europe entière, via Planetary Computer) — défaut
  # "ign" = RGE ALTI 1 m (Géoplateforme IGN, France métropolitaine uniquement)
  # "auto" = Copernicus par défaut, IGN si l'AOI est en France
  #
  # Note : les données TreeSatAI proviennent de Basse-Saxe (Allemagne),
  # donc seul Copernicus couvre la zone d'entraînement.
  # IGN est utile uniquement pour la prédiction spatiale sur une AOI française.
  dem_source = "copernicus",

  # URL WCS de la Géoplateforme IGN (MNT 1 m RGE ALTI, France uniquement)
  ign_wcs_url = "https://data.geopf.fr/wcs/ows",
  ign_coverage_id = "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES",

  # Features terrain à calculer
  compute_slope  = TRUE,
  compute_aspect = TRUE,
  compute_twi    = TRUE,
  compute_tpi    = TRUE,

  # Résolution de rééchantillonnage pour l'entraînement (mètres)
  resample_res = 10,

  # Plafond TWI (valeurs extrêmes dans les zones plates)
  twi_max = 20,

  # Rayon du filtre focal pour le TPI (en nombre de cellules)
  # TPI = DEM - DEM_lissé (moyenne focale)
  # 11 = fenêtre 11×11 pixels (110 m à 10 m de résolution)
  tpi_window = 11
)

# --- Masque forestier (OSO + NDVI) --------------------------------------------
# Masque combiné pour restreindre la classification aux zones boisées.
# Deux sources complémentaires :
#   1. OSO (CESBIO) : carte d'occupation du sol annuelle à 10 m (Sentinel-2)
#      Classes forestières : 31 = Forêt de feuillus, 32 = Forêt de conifères,
#                            33 = Forêt mixte (nomenclature OSO 2021+)
#   2. NDVI max annuel : seuil sur le maximum NDVI de la série temporelle
#      Détecte les coupes rases récentes (non encore visibles dans OSO)

FOREST_MASK_PARAMS <- list(
  # Activer le masque forestier combiné
  apply_forest_mask = TRUE,

  # --- OSO (carte d'occupation du sol CESBIO) ---
  use_oso = TRUE,

  # --- OCS GE (IGN, par département) ---
  # Téléchargement par département via API Géoplateforme (GPKG, quelques Mo)
  # API : https://data.geopf.fr/telechargement/resource/OCSGE
  # Ref : https://geoservices.ign.fr/ocsge#telechargement
  # Nomenclature couverture du sol (CS) :
  #   CS2.1.1.1 = Peuplements de feuillus
  #   CS2.1.1.2 = Peuplements de conifères
  #   CS2.1.1.3 = Peuplements mixtes
  ocsge_forest_pattern = "^CS2\\.1\\.1",

  # --- OSO raster CESBIO (France entière, ~6 Go) ---
  # Source : Recherche Data Gouv (CESBIO/CNES)
  # DOI : 10.57745/UZ2NJ7
  # Format : GeoTIFF, 10 m, Lambert-93 (EPSG:2154), France métropolitaine
  # Le raster France entière (~6 Go) est téléchargé dans un cache global
  # et découpé à l'AOI par projet (même logique que le package nemeton).
  oso_download_url = "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/8M1AN1",
  oso_manual_url   = "https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/UZ2NJ7",

  # Année (doit correspondre à l'année d'analyse)
  oso_year = 2021,

  # Classes OSO considérées comme forestières
  # Nomenclature OSO 23 classes (produit CESBIO/CNES depuis 2018) :
  #   Valeur pixel | Classe
  #   -------------|---------------------------
  #   1-4          | Urbain, routes
  #   5-12         | Cultures (colza, céréales, maïs, riz...)
  #   13           | Prairies
  #   14           | Vergers
  #   15           | Vignes
  #   16           | Forêts de feuillus  ← FORÊT
  #   17           | Forêts de conifères ← FORÊT
  #   18           | Pelouses naturelles
  #   19           | Landes ligneuses
  #   20           | Surfaces minérales / roche nue
  #   21           | Eau
  #   22           | Glaciers / neiges éternelles
  #   23           | Plages et dunes
  # Source : https://collections.sentinel-hub.com/cnes-land-cover-map/
  oso_forest_classes = c(16L, 17L),

  # Méthode de rééchantillonnage (catégoriel → nearest neighbor)
  oso_resample_method = "near",

  # --- NDVI max annuel ---
  use_ndvi = TRUE,

  # Seuil NDVI max : un pixel est considéré comme « végétation active » si

  # son NDVI max annuel dépasse ce seuil. Les zones en dessous (sol nu, urbain,
  # eau, coupes très récentes) sont exclues de la prédiction d'essence.
  # Valeur recommandée : 0.3–0.4
  ndvi_min_threshold = 0.4,

  # --- Combinaison des masques ---
  # "union"        : pixel boisé si OSO forêt OU NDVI > seuil
  #                  (conserve les coupes rases absentes d'OSO mais avec NDVI résiduel,
  #                   et les forêts récentes pas encore dans OSO)
  # "intersection" : pixel boisé si OSO forêt ET NDVI > seuil
  #                  (plus strict, exclut les coupes rases même si OSO dit forêt)
  combine_method = "union"
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
