#!/usr/bin/env python3
"""
TreeSatAI-Time-Series — Entraînement des modèles Deep Learning

Usage :
    conda activate treesat
    python python/train.py --data data/processed/feature_matrix.csv
    python python/train.py --data data/processed/feature_matrix.csv --model transformer --epochs 150
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR
from sklearn.metrics import (
    classification_report, confusion_matrix, accuracy_score, cohen_kappa_score
)

from config import (
    MODELS_DIR, OUTPUT_DIR, FIGURES_DIR,
    N_CLASSES, SPECIES_NAMES, TRAIN_CONFIG, N_S2_BANDS, TS_CONFIG
)
from models.architectures import get_model
from utils.dataset import create_dataloaders, TreeSatDataset


def train_one_epoch(model, loader, optimizer, criterion, device):
    """Entraîne le modèle sur une époque."""
    model.train()
    total_loss = 0
    correct = 0
    total = 0

    for batch in loader:
        x_s2 = batch["s2"].to(device)
        y = batch["label"].to(device)

        optimizer.zero_grad()

        if "s1" in batch:
            logits = model(x_s2, batch["s1"].to(device))
        else:
            logits = model(x_s2)

        loss = criterion(logits, y)
        loss.backward()

        # Gradient clipping
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)

        optimizer.step()

        total_loss += loss.item() * y.size(0)
        preds = logits.argmax(dim=1)
        correct += (preds == y).sum().item()
        total += y.size(0)

    return total_loss / total, correct / total


@torch.no_grad()
def evaluate(model, loader, criterion, device):
    """Évalue le modèle sur un dataset."""
    model.eval()
    total_loss = 0
    all_preds = []
    all_labels = []

    for batch in loader:
        x_s2 = batch["s2"].to(device)
        y = batch["label"].to(device)

        if "s1" in batch:
            logits = model(x_s2, batch["s1"].to(device))
        else:
            logits = model(x_s2)

        loss = criterion(logits, y)
        total_loss += loss.item() * y.size(0)

        all_preds.append(logits.argmax(dim=1).cpu())
        all_labels.append(y.cpu())

    all_preds = torch.cat(all_preds).numpy()
    all_labels = torch.cat(all_labels).numpy()
    avg_loss = total_loss / len(all_labels)
    accuracy = accuracy_score(all_labels, all_preds)

    return avg_loss, accuracy, all_preds, all_labels


def train(args):
    """Pipeline d'entraînement complet."""
    print("=" * 70)
    print("TreeSatAI-Time-Series — Entraînement Deep Learning")
    print("=" * 70)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device : {device}")
    if device.type == "cuda":
        print(f"  GPU : {torch.cuda.get_device_name(0)}")
        print(f"  VRAM : {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} Go")

    # --- 1. Données ---
    print(f"\n--- Chargement des données ---")
    print(f"Source : {args.data}")

    loaders = create_dataloaders(
        csv_path=args.data,
        batch_size=args.batch_size,
        train_ratio=TRAIN_CONFIG["train_ratio"],
        val_ratio=TRAIN_CONFIG["val_ratio"],
        seed=TRAIN_CONFIG["seed"],
        num_workers=args.workers
    )

    # Inférer les dimensions depuis le premier batch
    sample_batch = next(iter(loaders["train"]))
    n_channels = sample_batch["s2"].shape[1]
    n_timesteps = sample_batch["s2"].shape[2]
    has_s1 = "s1" in sample_batch

    print(f"Dimensions S2 : {n_channels} bandes × {n_timesteps} timesteps")
    if has_s1:
        print(f"Dimensions S1 : {sample_batch['s1'].shape[1]} bandes × {sample_batch['s1'].shape[2]} timesteps")

    # --- 2. Modèle ---
    print(f"\n--- Modèle : {args.model} ---")

    model_kwargs = {
        "n_channels": n_channels,
        "n_timesteps": n_timesteps,
        "n_classes": N_CLASSES,
        "dropout": args.dropout,
    }

    if args.model == "multisource" and has_s1:
        model_kwargs.update({
            "n_s2_channels": n_channels,
            "n_s2_timesteps": n_timesteps,
            "n_s1_channels": sample_batch["s1"].shape[1],
            "n_s1_timesteps": sample_batch["s1"].shape[2],
        })
        del model_kwargs["n_channels"]
        del model_kwargs["n_timesteps"]

    model = get_model(args.model, **model_kwargs).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Paramètres : {n_params:,}")
    print(model)

    # --- 3. Optimisation ---
    optimizer = AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = CosineAnnealingLR(optimizer, T_max=args.epochs, eta_min=args.lr * 0.01)

    # Poids des classes (pour données déséquilibrées)
    train_labels = []
    for batch in loaders["train"]:
        train_labels.append(batch["label"])
    train_labels = torch.cat(train_labels)
    class_counts = torch.bincount(train_labels, minlength=N_CLASSES).float()
    class_weights = (1.0 / (class_counts + 1)).to(device)
    class_weights = class_weights / class_weights.sum() * N_CLASSES

    criterion = nn.CrossEntropyLoss(weight=class_weights)

    # --- 4. Boucle d'entraînement ---
    print(f"\n--- Entraînement ({args.epochs} époques) ---")

    best_val_acc = 0
    patience_counter = 0
    history = {"train_loss": [], "train_acc": [], "val_loss": [], "val_acc": []}

    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    best_model_path = MODELS_DIR / f"treesatai_{args.model}_best.pt"

    t_start = time.time()

    for epoch in range(1, args.epochs + 1):
        # Train
        train_loss, train_acc = train_one_epoch(model, loaders["train"],
                                                  optimizer, criterion, device)
        # Validation
        val_loss, val_acc, _, _ = evaluate(model, loaders["val"], criterion, device)

        scheduler.step()

        history["train_loss"].append(train_loss)
        history["train_acc"].append(train_acc)
        history["val_loss"].append(val_loss)
        history["val_acc"].append(val_acc)

        # Early stopping
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            patience_counter = 0
            torch.save({
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "val_acc": val_acc,
                "model_name": args.model,
                "model_kwargs": model_kwargs,
                "class_names": SPECIES_NAMES,
            }, best_model_path)
        else:
            patience_counter += 1

        if epoch % 10 == 0 or epoch == 1:
            lr = optimizer.param_groups[0]["lr"]
            print(f"  Epoch {epoch:3d}/{args.epochs} | "
                  f"Train: loss={train_loss:.4f} acc={train_acc:.3f} | "
                  f"Val: loss={val_loss:.4f} acc={val_acc:.3f} | "
                  f"LR={lr:.6f} | Best={best_val_acc:.3f}")

        if patience_counter >= args.patience:
            print(f"\n  Early stopping à l'époque {epoch} (patience={args.patience})")
            break

    t_elapsed = time.time() - t_start
    print(f"\nEntraînement terminé en {t_elapsed/60:.1f} minutes")
    print(f"Meilleur modèle : {best_model_path} (val_acc={best_val_acc:.4f})")

    # --- 5. Évaluation sur le test set ---
    print(f"\n--- Évaluation finale (test set) ---")

    # Charger le meilleur modèle
    checkpoint = torch.load(best_model_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model_state_dict"])

    test_loss, test_acc, y_pred, y_true = evaluate(model, loaders["test"],
                                                     criterion, device)

    kappa = cohen_kappa_score(y_true, y_pred)
    report = classification_report(y_true, y_pred, target_names=SPECIES_NAMES,
                                    output_dict=True, zero_division=0)

    print(f"\n  Overall Accuracy : {test_acc:.4f} ({test_acc*100:.2f}%)")
    print(f"  Kappa de Cohen   : {kappa:.4f}")
    print(f"  Macro F1-Score   : {report['macro avg']['f1-score']:.4f}")

    print(f"\n  Classification par espèce :")
    print(classification_report(y_true, y_pred, target_names=SPECIES_NAMES,
                                 zero_division=0))

    # --- 6. Sauvegarder les résultats ---
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Historique
    import pandas as pd
    pd.DataFrame(history).to_csv(OUTPUT_DIR / f"{args.model}_training_history.csv", index=False)

    # Métriques par classe
    per_class = pd.DataFrame({
        "species": SPECIES_NAMES,
        "precision": [report[sp]["precision"] for sp in SPECIES_NAMES],
        "recall": [report[sp]["recall"] for sp in SPECIES_NAMES],
        "f1": [report[sp]["f1-score"] for sp in SPECIES_NAMES],
        "support": [report[sp]["support"] for sp in SPECIES_NAMES],
    })
    per_class.to_csv(OUTPUT_DIR / f"{args.model}_per_class_metrics.csv", index=False)

    # Matrice de confusion
    cm = confusion_matrix(y_true, y_pred)
    np.savetxt(OUTPUT_DIR / f"{args.model}_confusion_matrix.csv", cm,
               delimiter=",", fmt="%d")

    # Métadonnées JSON
    metadata = {
        "model": args.model,
        "n_params": n_params,
        "epochs_trained": epoch,
        "best_epoch": checkpoint["epoch"],
        "test_accuracy": float(test_acc),
        "kappa": float(kappa),
        "macro_f1": float(report["macro avg"]["f1-score"]),
        "training_time_min": round(t_elapsed / 60, 1),
        "device": str(device),
        "batch_size": args.batch_size,
        "learning_rate": args.lr,
    }
    with open(OUTPUT_DIR / f"{args.model}_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"\nRésultats sauvegardés dans {OUTPUT_DIR}")
    print(f"Modèle sauvegardé : {best_model_path}")

    return model, history, report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TreeSatAI-TS — Entraînement DL")

    parser.add_argument("--data", type=str, required=True,
                        help="Chemin vers feature_matrix.csv")
    parser.add_argument("--model", type=str, default="tempcnn",
                        choices=["tempcnn", "lstm", "transformer", "inception", "multisource"],
                        help="Architecture du modèle (défaut: tempcnn)")
    parser.add_argument("--epochs", type=int, default=TRAIN_CONFIG["n_epochs"],
                        help=f"Nombre d'époques (défaut: {TRAIN_CONFIG['n_epochs']})")
    parser.add_argument("--batch-size", type=int, default=TRAIN_CONFIG["batch_size"],
                        help=f"Taille de batch (défaut: {TRAIN_CONFIG['batch_size']})")
    parser.add_argument("--lr", type=float, default=TRAIN_CONFIG["learning_rate"],
                        help=f"Learning rate (défaut: {TRAIN_CONFIG['learning_rate']})")
    parser.add_argument("--weight-decay", type=float, default=TRAIN_CONFIG["weight_decay"])
    parser.add_argument("--dropout", type=float, default=0.3)
    parser.add_argument("--patience", type=int, default=TRAIN_CONFIG["patience"],
                        help=f"Early stopping patience (défaut: {TRAIN_CONFIG['patience']})")
    parser.add_argument("--workers", type=int, default=4,
                        help="Nombre de workers DataLoader")

    args = parser.parse_args()
    train(args)
