"""
TreeSatAI-Time-Series — Architectures Deep Learning (PyTorch)

Modèles pour la classification d'essences par séries temporelles :
  - TempCNN       : CNN 1D temporel (Pelletier et al., 2019)
  - LSTM_Classifier : LSTM bidirectionnel
  - TransformerTS : Transformer pour séries temporelles
  - InceptionTime : Inception adapté aux séries temporelles
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


# ==============================================================================
# 1. TempCNN — CNN temporel 1D (référence pour TreeSatAI)
# ==============================================================================

class TempCNN(nn.Module):
    """
    Temporal CNN pour classification de séries temporelles Sentinel-2.
    Architecture inspirée de Pelletier et al. (2019).

    Entrée : (batch, n_bands, n_timesteps)
    Sortie : (batch, n_classes) — logits

    3 blocs convolutifs avec BatchNorm + ReLU + Dropout,
    suivi d'un Global Average Pooling et d'une couche FC.
    """

    def __init__(self, n_channels=10, n_timesteps=73, n_classes=21,
                 filters=(64, 128, 256), kernel_sizes=(7, 5, 3),
                 dropout=0.3):
        super().__init__()

        self.n_channels = n_channels
        self.n_classes = n_classes

        # Bloc 1
        self.conv1 = nn.Conv1d(n_channels, filters[0], kernel_sizes[0],
                               padding=kernel_sizes[0] // 2)
        self.bn1 = nn.BatchNorm1d(filters[0])
        self.drop1 = nn.Dropout(dropout)

        # Bloc 2
        self.conv2 = nn.Conv1d(filters[0], filters[1], kernel_sizes[1],
                               padding=kernel_sizes[1] // 2)
        self.bn2 = nn.BatchNorm1d(filters[1])
        self.drop2 = nn.Dropout(dropout)

        # Bloc 3
        self.conv3 = nn.Conv1d(filters[1], filters[2], kernel_sizes[2],
                               padding=kernel_sizes[2] // 2)
        self.bn3 = nn.BatchNorm1d(filters[2])
        self.drop3 = nn.Dropout(dropout)

        # Classification head
        self.fc = nn.Linear(filters[2], n_classes)

    def forward(self, x):
        # x: (B, C, T)
        x = self.drop1(F.relu(self.bn1(self.conv1(x))))
        x = self.drop2(F.relu(self.bn2(self.conv2(x))))
        x = self.drop3(F.relu(self.bn3(self.conv3(x))))

        # Global Average Pooling sur la dimension temporelle
        x = x.mean(dim=2)  # (B, filters[-1])

        return self.fc(x)


# ==============================================================================
# 2. LSTM bidirectionnel
# ==============================================================================

class LSTMClassifier(nn.Module):
    """
    LSTM bidirectionnel pour séries temporelles.

    Entrée : (batch, n_bands, n_timesteps) → permutée en (batch, T, C)
    Sortie : (batch, n_classes)
    """

    def __init__(self, n_channels=10, n_timesteps=73, n_classes=21,
                 hidden_size=128, n_layers=2, dropout=0.3):
        super().__init__()

        self.lstm = nn.LSTM(
            input_size=n_channels,
            hidden_size=hidden_size,
            num_layers=n_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if n_layers > 1 else 0
        )

        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(hidden_size * 2, n_classes)  # *2 car bidirectionnel

    def forward(self, x):
        # x: (B, C, T) → (B, T, C)
        x = x.permute(0, 2, 1)

        output, (h_n, _) = self.lstm(x)

        # Concaténer les derniers hidden states des deux directions
        h_forward = h_n[-2]   # dernière couche, forward
        h_backward = h_n[-1]  # dernière couche, backward
        h_cat = torch.cat([h_forward, h_backward], dim=1)

        h_cat = self.dropout(h_cat)
        return self.fc(h_cat)


# ==============================================================================
# 3. Transformer pour séries temporelles
# ==============================================================================

class TransformerTS(nn.Module):
    """
    Transformer encoder pour classification de séries temporelles.

    Entrée : (batch, n_bands, n_timesteps)
    Sortie : (batch, n_classes)

    Positional encoding + TransformerEncoder + CLS token pooling.
    """

    def __init__(self, n_channels=10, n_timesteps=73, n_classes=21,
                 d_model=128, n_heads=4, n_layers=3, dropout=0.2):
        super().__init__()

        self.d_model = d_model

        # Projection des features spectrales vers d_model
        self.input_proj = nn.Linear(n_channels, d_model)

        # Positional encoding fixe (sinusoïdal)
        pe = torch.zeros(n_timesteps, d_model)
        position = torch.arange(0, n_timesteps, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-torch.log(torch.tensor(10000.0)) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe.unsqueeze(0))  # (1, T, d_model)

        # CLS token
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model))

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(d_model, n_classes)

    def forward(self, x):
        # x: (B, C, T) → (B, T, C)
        x = x.permute(0, 2, 1)
        B, T, C = x.shape

        # Projection
        x = self.input_proj(x)  # (B, T, d_model)

        # Ajouter positional encoding
        x = x + self.pe[:, :T, :]

        # Ajouter CLS token
        cls = self.cls_token.expand(B, -1, -1)
        x = torch.cat([cls, x], dim=1)  # (B, T+1, d_model)

        # Transformer
        x = self.transformer(x)

        # Prendre le CLS token
        cls_output = x[:, 0, :]

        cls_output = self.dropout(cls_output)
        return self.fc(cls_output)


# ==============================================================================
# 4. InceptionTime — Inception adapté aux séries temporelles
# ==============================================================================

class InceptionBlock(nn.Module):
    """Un bloc Inception pour séries temporelles."""

    def __init__(self, in_channels, n_filters=32, bottleneck_size=32):
        super().__init__()

        # Bottleneck
        self.bottleneck = nn.Conv1d(in_channels, bottleneck_size, 1)

        # 3 branches parallèles avec différentes tailles de kernel
        self.conv_10 = nn.Conv1d(bottleneck_size, n_filters, 10, padding=5)
        self.conv_20 = nn.Conv1d(bottleneck_size, n_filters, 20, padding=10)
        self.conv_40 = nn.Conv1d(bottleneck_size, n_filters, 40, padding=20)

        # Branche max pooling
        self.maxpool = nn.MaxPool1d(3, stride=1, padding=1)
        self.conv_pool = nn.Conv1d(in_channels, n_filters, 1)

        # BatchNorm sur la concaténation
        self.bn = nn.BatchNorm1d(n_filters * 4)

    def forward(self, x):
        bottleneck = self.bottleneck(x)

        branch1 = self.conv_10(bottleneck)
        branch2 = self.conv_20(bottleneck)
        branch3 = self.conv_40(bottleneck)
        branch4 = self.conv_pool(self.maxpool(x))

        # Aligner les tailles temporelles (au cas où)
        min_t = min(branch1.size(2), branch2.size(2), branch3.size(2), branch4.size(2))
        out = torch.cat([
            branch1[:, :, :min_t],
            branch2[:, :, :min_t],
            branch3[:, :, :min_t],
            branch4[:, :, :min_t]
        ], dim=1)

        return F.relu(self.bn(out))


class InceptionTime(nn.Module):
    """
    InceptionTime pour séries temporelles (Fawaz et al., 2020).

    Entrée : (batch, n_bands, n_timesteps)
    Sortie : (batch, n_classes)
    """

    def __init__(self, n_channels=10, n_timesteps=73, n_classes=21,
                 n_filters=32, n_blocks=3, dropout=0.2):
        super().__init__()

        self.blocks = nn.ModuleList()

        in_ch = n_channels
        for i in range(n_blocks):
            self.blocks.append(InceptionBlock(in_ch, n_filters))
            in_ch = n_filters * 4  # 4 branches concaténées

        self.gap = nn.AdaptiveAvgPool1d(1)
        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(n_filters * 4, n_classes)

    def forward(self, x):
        for block in self.blocks:
            x = block(x)

        x = self.gap(x).squeeze(-1)
        x = self.dropout(x)
        return self.fc(x)


# ==============================================================================
# 5. Multi-source : S2 + S1 (optionnel)
# ==============================================================================

class MultiSourceTempCNN(nn.Module):
    """
    TempCNN multi-source : branche S2 (optique) + branche S1 (radar).
    Fusion tardive par concaténation des features avant la couche FC.

    Entrée S2 : (batch, 10, 73)
    Entrée S1 : (batch, 3, 30)  — VV, VH, ratio × ~30 dates
    Sortie    : (batch, n_classes)
    """

    def __init__(self, n_s2_channels=10, n_s2_timesteps=73,
                 n_s1_channels=3, n_s1_timesteps=30,
                 n_classes=21, dropout=0.3):
        super().__init__()

        # Branche S2 (optique)
        self.s2_conv1 = nn.Conv1d(n_s2_channels, 64, 7, padding=3)
        self.s2_bn1 = nn.BatchNorm1d(64)
        self.s2_conv2 = nn.Conv1d(64, 128, 5, padding=2)
        self.s2_bn2 = nn.BatchNorm1d(128)
        self.s2_conv3 = nn.Conv1d(128, 128, 3, padding=1)
        self.s2_bn3 = nn.BatchNorm1d(128)

        # Branche S1 (radar)
        self.s1_conv1 = nn.Conv1d(n_s1_channels, 32, 5, padding=2)
        self.s1_bn1 = nn.BatchNorm1d(32)
        self.s1_conv2 = nn.Conv1d(32, 64, 3, padding=1)
        self.s1_bn2 = nn.BatchNorm1d(64)

        # Fusion
        self.dropout = nn.Dropout(dropout)
        self.fc1 = nn.Linear(128 + 64, 128)
        self.fc2 = nn.Linear(128, n_classes)

    def forward(self, x_s2, x_s1=None):
        # Branche S2
        s2 = F.relu(self.s2_bn1(self.s2_conv1(x_s2)))
        s2 = F.relu(self.s2_bn2(self.s2_conv2(s2)))
        s2 = F.relu(self.s2_bn3(self.s2_conv3(s2)))
        s2 = s2.mean(dim=2)  # GAP

        if x_s1 is not None:
            # Branche S1
            s1 = F.relu(self.s1_bn1(self.s1_conv1(x_s1)))
            s1 = F.relu(self.s1_bn2(self.s1_conv2(s1)))
            s1 = s1.mean(dim=2)  # GAP

            # Fusion tardive
            fused = torch.cat([s2, s1], dim=1)
        else:
            # S2 seul — padding pour la couche FC
            fused = torch.cat([s2, torch.zeros(s2.size(0), 64, device=s2.device)], dim=1)

        fused = self.dropout(fused)
        fused = F.relu(self.fc1(fused))
        fused = self.dropout(fused)
        return self.fc2(fused)


# ==============================================================================
# Factory
# ==============================================================================

def get_model(name="tempcnn", **kwargs):
    """Instanciation d'un modèle par nom."""
    models = {
        "tempcnn": TempCNN,
        "lstm": LSTMClassifier,
        "transformer": TransformerTS,
        "inception": InceptionTime,
        "multisource": MultiSourceTempCNN,
    }
    if name not in models:
        raise ValueError(f"Modèle inconnu '{name}'. Choix : {list(models.keys())}")
    return models[name](**kwargs)
