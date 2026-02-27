"""
TreeSatAI-Time-Series — Dataset et DataLoader PyTorch
"""

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from pathlib import Path
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import S2_BANDS, S1_BANDS, TS_CONFIG, SPECIES, N_CLASSES


class TreeSatDataset(Dataset):
    """
    Dataset PyTorch pour séries temporelles Sentinel-2 (+ S1 optionnel).

    Chaque échantillon = 1 parcelle avec :
      - X_s2 : (n_bands=10, n_timesteps=73) — séries temporelles S2
      - X_s1 : (n_bands=3, n_timesteps=30) — séries temporelles S1 (optionnel)
      - y    : int (0-19) — code espèce

    Peut charger depuis :
      - CSV (sortie du pipeline R)
      - Tenseur .pt
      - Matrice numpy
    """

    def __init__(self, data, labels=None, s1_data=None,
                 normalize=True, augment=False):
        """
        Args:
            data: np.ndarray (n_samples, n_bands, n_timesteps) ou chemin CSV
            labels: np.ndarray (n_samples,) — codes espèces (0-indexed)
            s1_data: np.ndarray (n_samples, 3, n_timesteps_s1) — S1 optionnel
            normalize: bool — normaliser les bandes
            augment: bool — augmentation temporelle aléatoire
        """
        if isinstance(data, (str, Path)):
            data, labels, s1_data = self._load_from_csv(data)

        self.X_s2 = torch.FloatTensor(data)
        self.y = torch.LongTensor(labels) if labels is not None else None
        self.X_s1 = torch.FloatTensor(s1_data) if s1_data is not None else None
        self.augment = augment
        self.has_s1 = self.X_s1 is not None

        if normalize:
            self._normalize()

    def _normalize(self):
        """Normalisation par bande (z-score sur l'ensemble du dataset)."""
        for b in range(self.X_s2.shape[1]):
            band = self.X_s2[:, b, :]
            mu = band.mean()
            sigma = band.std() + 1e-8
            self.X_s2[:, b, :] = (band - mu) / sigma

        if self.has_s1:
            for b in range(self.X_s1.shape[1]):
                band = self.X_s1[:, b, :]
                mu = band.mean()
                sigma = band.std() + 1e-8
                self.X_s1[:, b, :] = (band - mu) / sigma

    def _load_from_csv(self, csv_path):
        """Charge un CSV au format du pipeline R (feature_matrix.csv)."""
        df = pd.read_csv(csv_path)

        # Retrouver les colonnes de séries temporelles S2
        n_dates = TS_CONFIG["n_dates"]
        n_bands = len(S2_BANDS)

        # Les colonnes sont nommées : B02_d001, B02_d006, ...
        # Trouver les DOY labels
        doy_cols = sorted(set(
            col.split("_d")[1] for col in df.columns
            if "_d" in col and col.split("_d")[0] in S2_BANDS
        ))

        X_s2 = np.zeros((len(df), n_bands, len(doy_cols)), dtype=np.float32)
        for b, band in enumerate(S2_BANDS):
            for t, doy in enumerate(doy_cols):
                col = f"{band}_d{doy}"
                if col in df.columns:
                    X_s2[:, b, t] = df[col].fillna(0).values

        # Labels
        labels = None
        if "species_code" in df.columns:
            labels = df["species_code"].values.astype(int) - 1  # 0-indexed
        elif "species_name" in df.columns:
            name_to_idx = {SPECIES[i+1]["french"]: i for i in range(N_CLASSES)}
            labels = df["species_name"].map(name_to_idx).fillna(0).values.astype(int)

        # S1 (si les colonnes existent)
        s1_data = None
        s1_cols_vv = [c for c in df.columns if c.startswith("S1_VV_d")]
        if len(s1_cols_vv) > 0:
            s1_cols_vh = [c for c in df.columns if c.startswith("S1_VH_d")]
            s1_cols_ratio = [c for c in df.columns if c.startswith("S1_ratio_d")]
            n_s1_dates = len(s1_cols_vv)
            s1_data = np.zeros((len(df), 3, n_s1_dates), dtype=np.float32)
            for t, (cvv, cvh) in enumerate(zip(s1_cols_vv, s1_cols_vh)):
                s1_data[:, 0, t] = df[cvv].fillna(0).values
                s1_data[:, 1, t] = df[cvh].fillna(0).values
            if len(s1_cols_ratio) == n_s1_dates:
                for t, cr in enumerate(s1_cols_ratio):
                    s1_data[:, 2, t] = df[cr].fillna(0).values
            else:
                s1_data[:, 2, :] = s1_data[:, 0, :] - s1_data[:, 1, :]

        return X_s2, labels, s1_data

    def __len__(self):
        return len(self.X_s2)

    def __getitem__(self, idx):
        x_s2 = self.X_s2[idx]

        if self.augment:
            x_s2 = self._temporal_augment(x_s2)

        sample = {"s2": x_s2}

        if self.has_s1:
            sample["s1"] = self.X_s1[idx]

        if self.y is not None:
            sample["label"] = self.y[idx]

        return sample

    def _temporal_augment(self, x):
        """Augmentation temporelle : jitter + shift + scaling."""
        # Jitter (bruit gaussien)
        if torch.rand(1) < 0.5:
            x = x + torch.randn_like(x) * 0.02

        # Temporal shift (décalage de ±2 timesteps)
        if torch.rand(1) < 0.3:
            shift = torch.randint(-2, 3, (1,)).item()
            x = torch.roll(x, shifts=shift, dims=1)

        # Scaling (facteur aléatoire par bande)
        if torch.rand(1) < 0.3:
            scale = 0.9 + torch.rand(x.shape[0], 1) * 0.2  # [0.9, 1.1]
            x = x * scale

        return x


def create_dataloaders(csv_path=None, X_s2=None, y=None, X_s1=None,
                       batch_size=64, train_ratio=0.7, val_ratio=0.15,
                       seed=42, num_workers=4):
    """
    Création des DataLoaders train/val/test.

    Args:
        csv_path: chemin vers feature_matrix.csv (prioritaire)
        X_s2, y, X_s1: données numpy (si pas de CSV)
        batch_size: taille de batch
        train_ratio, val_ratio: proportions
        seed: graine aléatoire
        num_workers: workers pour le chargement

    Returns:
        dict avec 'train', 'val', 'test' DataLoaders
    """
    if csv_path is not None:
        full_dataset = TreeSatDataset(csv_path, normalize=False)
        X_s2_all = full_dataset.X_s2.numpy()
        y_all = full_dataset.y.numpy()
        X_s1_all = full_dataset.X_s1.numpy() if full_dataset.has_s1 else None
    else:
        X_s2_all = X_s2
        y_all = y
        X_s1_all = X_s1

    n_total = len(X_s2_all)
    np.random.seed(seed)

    # Split stratifié
    indices = np.arange(n_total)
    np.random.shuffle(indices)

    # Stratification par classe
    train_idx, val_idx, test_idx = [], [], []
    for cls in range(N_CLASSES):
        cls_idx = indices[y_all[indices] == cls]
        n_cls = len(cls_idx)
        n_train = int(n_cls * train_ratio)
        n_val = int(n_cls * val_ratio)

        train_idx.extend(cls_idx[:n_train])
        val_idx.extend(cls_idx[n_train:n_train + n_val])
        test_idx.extend(cls_idx[n_train + n_val:])

    # Créer les datasets
    def make_dataset(idx, augment=False):
        s1 = X_s1_all[idx] if X_s1_all is not None else None
        return TreeSatDataset(X_s2_all[idx], y_all[idx], s1,
                              normalize=True, augment=augment)

    train_ds = make_dataset(train_idx, augment=True)
    val_ds = make_dataset(val_idx, augment=False)
    test_ds = make_dataset(test_idx, augment=False)

    loaders = {
        "train": DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                            num_workers=num_workers, pin_memory=True),
        "val":   DataLoader(val_ds, batch_size=batch_size, shuffle=False,
                            num_workers=num_workers, pin_memory=True),
        "test":  DataLoader(test_ds, batch_size=batch_size, shuffle=False,
                            num_workers=num_workers, pin_memory=True),
    }

    print(f"Datasets : train={len(train_ds)}, val={len(val_ds)}, test={len(test_ds)}")
    return loaders
