#!/usr/bin/env python3
"""
TreeSatAI-Time-Series — Prédiction pixel par pixel avec modèle PyTorch

Usage depuis R (via reticulate) ou en ligne de commande :
    python python/predict.py --model output/models/treesatai_tempcnn_best.pt \
                             --input data/processed/pixel_features.npy \
                             --output output/predictions.npy
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import N_CLASSES, SPECIES_NAMES, MODELS_DIR
from models.architectures import get_model


def load_model(model_path, device=None):
    """
    Charge un modèle PyTorch sauvegardé.

    Args:
        model_path: chemin vers le fichier .pt
        device: 'cuda' ou 'cpu'

    Returns:
        model, class_names, metadata
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    checkpoint = torch.load(model_path, map_location=device, weights_only=False)

    model_name = checkpoint.get("model_name", "tempcnn")
    model_kwargs = checkpoint.get("model_kwargs", {})
    class_names = checkpoint.get("class_names", SPECIES_NAMES)

    model = get_model(model_name, **model_kwargs)
    model.load_state_dict(checkpoint["model_state_dict"])
    model = model.to(device)
    model.eval()

    print(f"Modèle chargé : {model_name} ({sum(p.numel() for p in model.parameters()):,} params)")
    print(f"  Val accuracy : {checkpoint.get('val_acc', 'N/A')}")
    print(f"  Classes : {len(class_names)}")

    return model, class_names, checkpoint


@torch.no_grad()
def predict_pixels(model, pixel_data, batch_size=512, device=None):
    """
    Prédiction sur un ensemble de pixels.

    Args:
        model: modèle PyTorch en mode eval
        pixel_data: np.ndarray (n_pixels, n_bands, n_timesteps)
        batch_size: taille de batch pour l'inférence
        device: device PyTorch

    Returns:
        dict avec:
          - predicted_class: np.ndarray (n_pixels,) — index de classe (0-based)
          - probabilities: np.ndarray (n_pixels, n_classes) — probabilités softmax
          - max_proba: np.ndarray (n_pixels,) — confiance max
    """
    if device is None:
        device = next(model.parameters()).device

    n_pixels = len(pixel_data)
    all_probs = []

    for start in range(0, n_pixels, batch_size):
        end = min(start + batch_size, n_pixels)
        batch = torch.FloatTensor(pixel_data[start:end]).to(device)

        logits = model(batch)
        probs = torch.softmax(logits, dim=1)
        all_probs.append(probs.cpu().numpy())

    all_probs = np.concatenate(all_probs, axis=0)
    predicted_class = all_probs.argmax(axis=1)
    max_proba = all_probs.max(axis=1)

    return {
        "predicted_class": predicted_class,
        "probabilities": all_probs,
        "max_proba": max_proba,
    }


@torch.no_grad()
def predict_multisource(model, s2_data, s1_data=None, batch_size=512, device=None):
    """
    Prédiction multi-source (S2 + S1).

    Args:
        model: modèle MultiSourceTempCNN
        s2_data: np.ndarray (n_pixels, n_s2_bands, n_s2_timesteps)
        s1_data: np.ndarray (n_pixels, n_s1_bands, n_s1_timesteps) ou None
        batch_size: taille de batch
        device: device PyTorch

    Returns:
        dict comme predict_pixels
    """
    if device is None:
        device = next(model.parameters()).device

    n_pixels = len(s2_data)
    all_probs = []

    for start in range(0, n_pixels, batch_size):
        end = min(start + batch_size, n_pixels)
        batch_s2 = torch.FloatTensor(s2_data[start:end]).to(device)
        batch_s1 = None
        if s1_data is not None:
            batch_s1 = torch.FloatTensor(s1_data[start:end]).to(device)

        logits = model(batch_s2, batch_s1)
        probs = torch.softmax(logits, dim=1)
        all_probs.append(probs.cpu().numpy())

    all_probs = np.concatenate(all_probs, axis=0)

    return {
        "predicted_class": all_probs.argmax(axis=1),
        "probabilities": all_probs,
        "max_proba": all_probs.max(axis=1),
    }


def export_onnx(model, output_path, n_channels=10, n_timesteps=73):
    """
    Export du modèle en ONNX pour utilisation hors PyTorch.

    Args:
        model: modèle PyTorch
        output_path: chemin de sortie .onnx
        n_channels: nombre de bandes
        n_timesteps: nombre de timesteps
    """
    model.eval()
    dummy_input = torch.randn(1, n_channels, n_timesteps)

    torch.onnx.export(
        model,
        dummy_input,
        str(output_path),
        input_names=["time_series"],
        output_names=["logits"],
        dynamic_axes={
            "time_series": {0: "batch_size"},
            "logits": {0: "batch_size"},
        },
        opset_version=17,
    )
    print(f"Modèle ONNX exporté : {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TreeSatAI-TS — Prédiction DL")

    parser.add_argument("--model", type=str, required=True,
                        help="Chemin vers le modèle .pt")
    parser.add_argument("--input", type=str, required=True,
                        help="Données d'entrée (.npy : n_pixels × n_bands × n_timesteps)")
    parser.add_argument("--output", type=str, default="output/predictions.npy",
                        help="Fichier de sortie")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--export-onnx", action="store_true",
                        help="Exporter aussi en ONNX")

    args = parser.parse_args()

    # Charger le modèle
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, class_names, meta = load_model(args.model, device)

    # Charger les données
    print(f"Chargement : {args.input}")
    pixel_data = np.load(args.input)
    print(f"  Shape : {pixel_data.shape}")

    # Prédiction
    results = predict_pixels(model, pixel_data, batch_size=args.batch_size, device=device)

    # Sauvegarder
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    np.save(output_path, results["predicted_class"])
    np.save(output_path.with_suffix(".proba.npy"), results["probabilities"])
    np.save(output_path.with_suffix(".confidence.npy"), results["max_proba"])

    print(f"\nPrédictions sauvegardées :")
    print(f"  Classes    : {output_path}")
    print(f"  Probas     : {output_path.with_suffix('.proba.npy')}")
    print(f"  Confiance  : {output_path.with_suffix('.confidence.npy')}")

    # Résumé
    unique, counts = np.unique(results["predicted_class"], return_counts=True)
    print(f"\nEspèces détectées :")
    for cls, cnt in sorted(zip(unique, counts), key=lambda x: -x[1]):
        pct = cnt / len(results["predicted_class"]) * 100
        name = class_names[cls] if cls < len(class_names) else f"Classe {cls}"
        print(f"  {name:30s} : {cnt:6d} pixels ({pct:5.1f}%)")

    # Export ONNX
    if args.export_onnx:
        onnx_path = Path(args.model).with_suffix(".onnx")
        export_onnx(model.cpu(), onnx_path,
                     n_channels=pixel_data.shape[1],
                     n_timesteps=pixel_data.shape[2])
