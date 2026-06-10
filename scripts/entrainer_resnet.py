"""
Entraînement d'un ResNet-18 sur les patches Landsat 7 pour prédire
la richesse moyenne d'un cluster DHS.

Usage:
    python scripts/entrainer_resnet.py

Ce script est conçu pour tourner sur Kaggle (GPU disponible).
Paquets requis : torch, torchvision, numpy, pandas, scikit-learn, matplotlib
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

PATCH_DIR = "data/processed/patches_landsat"
LABELS_CSV = os.path.join(PATCH_DIR, "labels.csv")
OUT_METRICS = "outputs/tables/cnn_metrics.csv"
OUT_PREDS = "outputs/tables/cnn_predictions.csv"
OUT_MODEL = "outputs/models/resnet18.pt"
OUT_FIG = "outputs/figures/cnn_pred_vs_true.png"

BATCH_SIZE = 32
LR = 5e-4
WEIGHT_DECAY = 1e-4
PATIENCE = 10
MAX_EPOCHS = 100
TEST_SIZE = 0.2


def main():
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
    import torchvision.models as models
    import torchvision.transforms as T
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import r2_score, mean_squared_error

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device : {device}")

    os.makedirs("outputs/tables", exist_ok=True)
    os.makedirs("outputs/models", exist_ok=True)
    os.makedirs("outputs/figures", exist_ok=True)

    # 1. Charger les labels
    labels = pd.read_csv(LABELS_CSV)
    labels["cluster"] = labels["cluster"].astype(int)
    print(f"Labels chargés : {len(labels)}")

    # Vérifier quels patches existent
    labels["patch_exists"] = labels["cluster"].apply(
        lambda c: os.path.exists(os.path.join(PATCH_DIR, f"{c}.npy"))
    )
    labels = labels[labels["patch_exists"]].reset_index(drop=True)
    print(f"Clusters avec patch : {len(labels)}")

    if len(labels) < 10:
        sys.exit("Pas assez de patches. Lancer d'abord telecharger_landsat.py")

    # 2. Dataset PyTorch
    class LandsatDataset(Dataset):
        def __init__(self, df, patch_dir, transform=None):
            self.df = df.reset_index(drop=True)
            self.patch_dir = patch_dir
            self.transform = transform

        def __len__(self):
            return len(self.df)

        def __getitem__(self, idx):
            row = self.df.iloc[idx]
            patch = np.load(os.path.join(
                self.patch_dir, f"{int(row['cluster'])}.npy"
            ))
            # (H, W, C) -> (C, H, W)
            patch = torch.from_numpy(patch).permute(2, 0, 1).float()
            if self.transform:
                patch = self.transform(patch)
            return patch, torch.tensor(row["wealth_mean"], dtype=torch.float32)

    # 3. Split train/test
    train_df, test_df = train_test_split(
        labels, test_size=TEST_SIZE, random_state=42
    )
    print(f"Train : {len(train_df)} | Test : {len(test_df)}")

    # 4. Transformations
    # Normalisation ImageNet (moyenne/écart-type des 3 canaux RGB)
    IMAGENET_MEAN = [0.485, 0.456, 0.406]
    IMAGENET_STD = [0.229, 0.224, 0.225]

    train_transform = T.Compose([
        T.Normalize(IMAGENET_MEAN, IMAGENET_STD),
        T.RandomHorizontalFlip(p=0.5),
        T.RandomRotation(15),
    ])
    test_transform = T.Compose([
        T.Normalize(IMAGENET_MEAN, IMAGENET_STD),
    ])

    train_ds = LandsatDataset(train_df, PATCH_DIR, transform=train_transform)
    test_ds = LandsatDataset(test_df, PATCH_DIR, transform=test_transform)

    train_loader = DataLoader(
        train_ds, batch_size=BATCH_SIZE, shuffle=True, num_workers=2
    )
    test_loader = DataLoader(
        test_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=2
    )

    # 5. Modèle ResNet-18
    model = models.resnet18(weights="DEFAULT")

    # Geler les premières couches
    for name, param in model.named_parameters():
        if "layer4" in name or "layer3" in name or "fc" in name:
            param.requires_grad = True
        else:
            param.requires_grad = False

    # Remplacer la dernière couche
    in_features = model.fc.in_features
    model.fc = nn.Sequential(
        nn.Linear(in_features, 256),
        nn.ReLU(),
        nn.Dropout(0.3),
        nn.Linear(256, 1),
    )
    model = model.to(device)

    # 6. Optimiseur (ne toucher que les paramètres entraînables)
    trainable_params = [p for p in model.parameters() if p.requires_grad]
    optimizer = optim.Adam(trainable_params, lr=LR, weight_decay=WEIGHT_DECAY)
    criterion = nn.MSELoss()
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=5
    )

    # 7. Boucle d'entraînement
    best_loss = float("inf")
    patience_counter = 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        train_loss = 0.0
        for x, y in train_loader:
            x, y = x.to(device), y.to(device)
            optimizer.zero_grad()
            pred = model(x).squeeze(-1)
            loss = criterion(pred, y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * x.size(0)
        train_loss /= len(train_ds)

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for x, y in test_loader:
                x, y = x.to(device), y.to(device)
                pred = model(x).squeeze(-1)
                loss = criterion(pred, y)
                val_loss += loss.item() * x.size(0)
        val_loss /= len(test_ds)

        scheduler.step(val_loss)

        if epoch % 5 == 0 or epoch == 1:
            print(f"Epoch {epoch:3d} | train MSE = {train_loss:.4f} | val MSE = {val_loss:.4f}")

        # Early stopping
        if val_loss < best_loss:
            best_loss = val_loss
            torch.save(model.state_dict(), OUT_MODEL)
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= PATIENCE:
                print(f"Early stopping à l'epoch {epoch}")
                break

    # 8. Évaluation finale
    model.load_state_dict(torch.load(OUT_MODEL))
    model.eval()
    all_preds, all_true = [], []
    with torch.no_grad():
        for x, y in test_loader:
            x = x.to(device)
            pred = model(x).squeeze(-1).cpu().numpy()
            all_preds.extend(pred)
            all_true.extend(y.numpy())

    all_preds = np.array(all_preds)
    all_true = np.array(all_true)

    rmse = float(np.sqrt(mean_squared_error(all_true, all_preds)))
    r2 = float(r2_score(all_true, all_preds))
    print(f"\nTest RMSE : {rmse:.3f}")
    print(f"Test R²   : {r2:.3f}")

    # 9. Sauvegardes
    pd.DataFrame([
        {"metric": "rmse_test", "value": rmse},
        {"metric": "r2_test", "value": r2},
        {"metric": "n_train", "value": len(train_df)},
        {"metric": "n_test", "value": len(test_df)},
    ]).to_csv(OUT_METRICS, index=False)

    pd.DataFrame({
        "cluster": test_df["cluster"].values,
        "wealth_true": all_true,
        "wealth_pred": all_preds,
    }).to_csv(OUT_PREDS, index=False)

    # 10. Figure
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(all_true, all_preds, alpha=0.7, edgecolor="k")
    lo = min(all_true.min(), all_preds.min())
    hi = max(all_true.max(), all_preds.max())
    ax.plot([lo, hi], [lo, hi], "r--", linewidth=1)
    ax.set_xlabel("Richesse vraie (DHS)")
    ax.set_ylabel("Richesse prédite (ResNet-18)")
    ax.set_title(f"ResNet-18 sur Landsat 7 : RMSE = {rmse:.1f}, R² = {r2:.3f}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=200)
    plt.close(fig)

    print(f"Métriques : {OUT_METRICS}")
    print(f"Prédictions : {OUT_PREDS}")
    print(f"Modèle : {OUT_MODEL}")
    print(f"Figure : {OUT_FIG}")

    # 11. Prédictions sur TOUS les clusters (pour le FH)
    full_ds = LandsatDataset(labels, PATCH_DIR, transform=test_transform)
    full_loader = DataLoader(
        full_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=2
    )
    all_full_preds = []
    model.eval()
    with torch.no_grad():
        for x, _ in full_loader:
            x = x.to(device)
            pred = model(x).squeeze(-1).cpu().numpy()
            all_full_preds.extend(pred)

    pd.DataFrame({
        "cluster": labels["cluster"].values,
        "cnn_prediction": all_full_preds,
    }).to_csv("data/processed/cnn_predictions_cluster.csv", index=False)
    print("Prédictions tous clusters : data/processed/cnn_predictions_cluster.csv")


if __name__ == "__main__":
    main()
