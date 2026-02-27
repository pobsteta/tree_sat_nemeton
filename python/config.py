"""
TreeSatAI-Time-Series — Configuration Python
"""

from pathlib import Path

# Chemins du projet
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
TS_DIR = DATA_DIR / "timeseries"
OUTPUT_DIR = PROJECT_ROOT / "output"
MODELS_DIR = OUTPUT_DIR / "models"
FIGURES_DIR = PROJECT_ROOT / "figures"

# Créer les dossiers
for d in [DATA_DIR, RAW_DIR, PROCESSED_DIR, TS_DIR, OUTPUT_DIR, MODELS_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# --- Les 20 espèces européennes -----------------------------------------------
SPECIES = {
    1:  {"latin": "Quercus robur",          "french": "Chêne pédonculé",      "type": "feuillu",  "phenology": "deciduous"},
    2:  {"latin": "Quercus petraea",         "french": "Chêne sessile",        "type": "feuillu",  "phenology": "deciduous"},
    3:  {"latin": "Quercus pubescens",       "french": "Chêne pubescent",      "type": "feuillu",  "phenology": "deciduous"},
    4:  {"latin": "Quercus ilex",            "french": "Chêne vert",           "type": "feuillu",  "phenology": "evergreen"},
    5:  {"latin": "Fagus sylvatica",         "french": "Hêtre",                "type": "feuillu",  "phenology": "deciduous"},
    6:  {"latin": "Castanea sativa",         "french": "Châtaignier",          "type": "feuillu",  "phenology": "deciduous"},
    7:  {"latin": "Carpinus betulus",        "french": "Charme",               "type": "feuillu",  "phenology": "deciduous"},
    8:  {"latin": "Betula pendula",          "french": "Bouleau verruqueux",   "type": "feuillu",  "phenology": "deciduous"},
    9:  {"latin": "Fraxinus excelsior",      "french": "Frêne commun",         "type": "feuillu",  "phenology": "deciduous"},
    10: {"latin": "Acer pseudoplatanus",     "french": "Érable sycomore",      "type": "feuillu",  "phenology": "deciduous"},
    11: {"latin": "Populus spp.",            "french": "Peupliers",            "type": "feuillu",  "phenology": "deciduous"},
    12: {"latin": "Robinia pseudoacacia",    "french": "Robinier faux-acacia", "type": "feuillu",  "phenology": "deciduous"},
    13: {"latin": "Picea abies",             "french": "Épicéa commun",        "type": "résineux", "phenology": "evergreen"},
    14: {"latin": "Abies alba",              "french": "Sapin pectiné",        "type": "résineux", "phenology": "evergreen"},
    15: {"latin": "Pseudotsuga menziesii",   "french": "Douglas",              "type": "résineux", "phenology": "evergreen"},
    16: {"latin": "Pinus sylvestris",        "french": "Pin sylvestre",        "type": "résineux", "phenology": "evergreen"},
    17: {"latin": "Pinus pinaster",          "french": "Pin maritime",         "type": "résineux", "phenology": "evergreen"},
    18: {"latin": "Pinus nigra",             "french": "Pin noir",             "type": "résineux", "phenology": "evergreen"},
    19: {"latin": "Pinus halepensis",        "french": "Pin d'Alep",           "type": "résineux", "phenology": "evergreen"},
    20: {"latin": "Larix decidua",           "french": "Mélèze d'Europe",      "type": "résineux", "phenology": "deciduous"},
}

SPECIES_NAMES = [SPECIES[i]["french"] for i in range(1, 21)]
N_CLASSES = len(SPECIES)

# Bandes Sentinel-2
S2_BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
N_S2_BANDS = len(S2_BANDS)

# Bandes Sentinel-1
S1_BANDS = ["VV", "VH"]
N_S1_BANDS = len(S1_BANDS)

# Paramètres séries temporelles
TS_CONFIG = {
    "target_interval_days": 5,
    "n_dates": 73,  # 365 / 5
    "sg_window": 7,
    "sg_order": 3,
}

# Paramètres d'entraînement par défaut
TRAIN_CONFIG = {
    "batch_size": 64,
    "learning_rate": 1e-3,
    "weight_decay": 1e-4,
    "n_epochs": 100,
    "patience": 15,       # early stopping
    "train_ratio": 0.7,
    "val_ratio": 0.15,
    "seed": 42,
}
