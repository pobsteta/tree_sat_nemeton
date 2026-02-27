# TreeSatAI-Time-Series (IGNF) — Classification d'essences forestières

Classification de **20 espèces forestières européennes** par analyse de
**séries temporelles Sentinel-2 annuelles**, basée sur les signatures
phénologiques (verdissement printanier, plateau estival, sénescence automnale).

## Principe

Chaque espèce d'arbre possède une **signature phénologique** unique : la façon
dont elle verdit au printemps et rougit à l'automne. En capturant cette dynamique
à travers une année complète d'images Sentinel-2 (tous les 5 jours), on peut
distinguer les espèces par leur "empreinte temporelle".

```
NDVI ▲
 0.8 │          ╭──── Hêtre (pic tardif, forte amplitude)
     │        ╭╯  ╲
 0.6 │  ═══════════════ Épicéa (stable, sempervirent)
     │      ╱        ╲
 0.4 │    ╱            ╲ ← Mélèze (résineux caduc !)
     │  ╱                ╲
 0.2 │╱                    ╲
     └──────────────────────────► Mois
      J  F  M  A  M  J  J  A  S  O  N  D
```

## Les 20 espèces

### Feuillus (12 espèces)
| # | Nom latin | Nom français | Phénologie |
|---|-----------|-------------|------------|
| 1 | *Quercus robur* | Chêne pédonculé | Caducifolié |
| 2 | *Quercus petraea* | Chêne sessile | Caducifolié |
| 3 | *Quercus pubescens* | Chêne pubescent | Caducifolié |
| 4 | *Quercus ilex* | Chêne vert | **Sempervirent** |
| 5 | *Fagus sylvatica* | Hêtre | Caducifolié |
| 6 | *Castanea sativa* | Châtaignier | Caducifolié |
| 7 | *Carpinus betulus* | Charme | Caducifolié |
| 8 | *Betula pendula* | Bouleau verruqueux | Caducifolié |
| 9 | *Fraxinus excelsior* | Frêne commun | Caducifolié |
| 10 | *Acer pseudoplatanus* | Érable sycomore | Caducifolié |
| 11 | *Populus spp.* | Peupliers | Caducifolié |
| 12 | *Robinia pseudoacacia* | Robinier faux-acacia | Caducifolié |

### Résineux (8 espèces)
| # | Nom latin | Nom français | Phénologie |
|---|-----------|-------------|------------|
| 13 | *Picea abies* | Épicéa commun | Sempervirent |
| 14 | *Abies alba* | Sapin pectiné | Sempervirent |
| 15 | *Pseudotsuga menziesii* | Douglas | Sempervirent |
| 16 | *Pinus sylvestris* | Pin sylvestre | Sempervirent |
| 17 | *Pinus pinaster* | Pin maritime | Sempervirent |
| 18 | *Pinus nigra* | Pin noir | Sempervirent |
| 19 | *Pinus halepensis* | Pin d'Alep | Sempervirent |
| 20 | *Larix decidua* | Mélèze d'Europe | **Caducifolié** |

> Le **Chêne vert** est le seul feuillu sempervirent du jeu de données.
> Le **Mélèze** est le seul résineux caducifolié (il perd ses aiguilles en hiver).

## Structure du projet

```
tree_sat_nemeton/
├── R/
│   ├── 00_config.R            # Configuration, paramètres, espèces
│   ├── 01_utils.R             # Fonctions utilitaires (indices, interpolation)
│   ├── 02_phenology.R         # Extraction de métriques phénologiques
│   ├── 03_data_acquisition.R  # Acquisition données (IFN, Sentinel-2, STAC)
│   ├── 04_visualization.R     # Profils, heatmaps, matrices de confusion
│   ├── 05_classification.R    # Random Forest + CNN temporel
│   └── 06_pipeline.R          # Pipeline complet
├── data/
│   ├── raw/                   # Données brutes
│   ├── processed/             # Matrice de features
│   └── timeseries/            # Séries temporelles extraites
├── output/
│   └── models/                # Modèles entraînés
├── figures/                   # Graphiques générés
└── README.md
```

## Utilisation

### Mode démonstration (données synthétiques)

```bash
cd tree_sat_nemeton
Rscript R/06_pipeline.R
```

Génère des profils phénologiques réalistes pour les 20 espèces, entraîne un
Random Forest et produit les visualisations.

### Mode données réelles

```bash
Rscript R/06_pipeline.R --data /chemin/vers/dataset
```

### Options

| Argument | Description | Défaut |
|----------|-------------|--------|
| `--data <path>` | Chemin vers le dataset TreeSatAI-TS | Synthétique |
| `--mode rf\|cnn\|both` | Méthode de classification | `rf` |
| `--samples <n>` | Échantillons par espèce (mode synthétique) | `50` |
| `--no-viz` | Désactiver les visualisations | `FALSE` |

### Avec CNN temporel (nécessite torch)

```bash
# Installer torch pour R
Rscript -e "install.packages('torch'); torch::install_torch()"

# Lancer avec le CNN
Rscript R/06_pipeline.R --mode both
```

## Features extraites

Le pipeline extrait **4 familles de features** à partir des séries temporelles :

### 1. Séries temporelles brutes (730 features)
- 10 bandes Sentinel-2 × 73 dates = 730 valeurs interpolées à 5 jours

### 2. Indices spectraux (6 indices × 73 dates + statistiques)
- **NDVI** — Normalized Difference Vegetation Index
- **EVI** — Enhanced Vegetation Index
- **NDWI** — Normalized Difference Water Index
- **CRI** — Carotenoid Reflectance Index
- **RENDVI** — Red Edge NDVI
- **NBR** — Normalized Burn Ratio

### 3. Métriques phénologiques (26 features)
- **SOS** (Start of Season) — début du verdissement
- **EOS** (End of Season) — fin de la sénescence
- **LOS** (Length of Season) — durée de la saison de croissance
- **Amplitude** — différence NDVI max - min
- **Taux de verdissement** — pente printanière
- **Taux de sénescence** — pente automnale
- **Intégrale saisonnière** — proxy de la productivité
- **NDVI saisonnier** — moyennes hiver/printemps/été/automne
- **Ratio été/hiver** — indicateur persistant vs caduc
- **Asymétrie** — rapport verdissement/sénescence

### 4. Harmoniques de Fourier (7 features)
- Amplitude et phase des 3 premières harmoniques
- Captent la périodicité et la forme du signal annuel

## Méthodes de classification

### Random Forest (ranger)
- 500 arbres, mtry = sqrt(n_features)
- Importance par permutation
- Validation croisée 5-fold

### CNN temporel 1D (torch, optionnel)
- Architecture : Conv1D(64) → Conv1D(128) → Conv1D(256) → GAP → FC(128) → FC(20)
- Entrée : 10 bandes × 73 pas de temps
- BatchNorm + Dropout (0.3)
- Optimiseur Adam, 50 époques

## Visualisations générées

| Fichier | Description |
|---------|-------------|
| `01_phenological_profiles_NDVI.png` | Profils NDVI moyens des 20 espèces |
| `02_deciduous_vs_evergreen.png` | Comparaison caducifolié vs sempervirent |
| `03_profiles_by_type.png` | Profils par type (feuillu/résineux) et phénologie |
| `04_confusion_matrix.png` | Matrice de confusion normalisée |
| `05_variable_importance.png` | Top 30 features les plus discriminantes |
| `06_species_metrics.png` | Précision / Rappel / F1 par espèce |
| `07_spectral_heatmap.png` | Heatmap des signatures NDVI (z-score) |
| `08_dashboard.png` | Tableau de bord récapitulatif |

## Données

### Source : IGN France / TreeSatAI

Le dataset TreeSatAI-Time-Series est construit à partir de :
- **Placettes IFN** (Inventaire Forestier National) : localisation et identification des espèces
- **Sentinel-2 L2A** : images satellite tous les 5 jours, 10 bandes spectrales
- **Période** : 1 année complète (janvier à décembre)

### Preprocessing
1. Masquage nuageux (SCL : classes 4, 5, 6)
2. Interpolation temporelle (spline naturelle → grille régulière à 5 jours)
3. Lissage Savitzky-Golay (polynôme ordre 3, fenêtre 7)
4. Normalisation et extraction de features

## Références

- Ahlswede, S. et al. (2023). "TreeSatAI Benchmark Archive". *Earth System Science Data*.
- Inglada, J. et al. (2017). "Operational High Resolution Land Cover Map Production at the Country Scale Using Satellite Image Time Series". *Remote Sensing*.
- IGN France — Inventaire Forestier National (IFN)
- Copernicus Sentinel-2 Data Space Ecosystem
