"""
Entraine un ResNet-18 sur les patches Landsat 7 pour predire
la richesse moyenne d'un cluster DHS (Niger).

Validation : Group 5-fold par departement (test spatial).
Courbes d'apprentissage + boxplot des R².

Usage:
    python scripts/entrainer_resnet.py
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
CLUSTERS_CSV = "data/processed/clusters_gps.csv"
OUT_METRICS = "outputs/tables/cnn_metrics.csv"
OUT_PREDS = "outputs/tables/cnn_predictions.csv"
OUT_MODEL = "outputs/models/resnet18.pt"
OUT_FIG_PRED = "outputs/figures/cnn_pred_vs_true.png"
OUT_FIG_CURVES = "outputs/figures/cnn_learning_curves.png"
OUT_FIG_CV = "outputs/figures/cnn_cv_boxplot.png"
OUT_HISTORY = "outputs/tables/cnn_training_history.csv"

BATCH_SIZE = 32
LR = 5e-4
WEIGHT_DECAY = 1e-4
PATIENCE = 10
MAX_EPOCHS = 100
N_FOLDS = 5


def load_patch(cluster_id, patch_dir):
    # Cherche d'abord .npy (deja normalise), puis .tif (brut)
    npy_path = os.path.join(patch_dir, f"{int(cluster_id)}.npy")
    tif_path = os.path.join(patch_dir, f"cluster_{int(cluster_id):03d}.tif")

    if os.path.exists(npy_path):
        return np.load(npy_path).astype(np.float32)

    if os.path.exists(tif_path):
        import rasterio
        with rasterio.open(tif_path) as src:
            patch = src.read()
        # (C, H, W) -> (H, W, C)
        patch = np.transpose(patch, (1, 2, 0)).astype(np.float32)
        # Normalisation percentile [0, 1]
        p2, p98 = np.percentile(patch, [2, 98])
        patch = np.clip((patch - p2) / (p98 - p2 + 1e-8), 0, 1)
        return patch

    return None


def train_one_fold(train_idx, test_idx, labels_df, patch_dir, device, fold):
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
    import torchvision.models as models
    import torchvision.transforms as T

    train_df = labels_df.iloc[train_idx].reset_index(drop=True)
    test_df = labels_df.iloc[test_idx].reset_index(drop=True)

    class LandsatDataset(Dataset):
        def __init__(self, df, patch_dir, transform=None):
            self.df = df
            self.patch_dir = patch_dir
            self.transform = transform

        def __len__(self):
            return len(self.df)

        def __getitem__(self, idx):
            row = self.df.iloc[idx]
            patch = load_patch(row["cluster"], self.patch_dir)
            if patch is None:
                raise FileNotFoundError(f"Patch {row['cluster']} introuvable")
            # (H, W, C) -> (C, H, W)
            patch = torch.from_numpy(patch).permute(2, 0, 1).float()
            if self.transform:
                patch = self.transform(patch)
            return patch, torch.tensor(row["wealth_mean"], dtype=torch.float32)

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

    train_ds = LandsatDataset(train_df, patch_dir, transform=train_transform)
    test_ds = LandsatDataset(test_df, patch_dir, transform=test_transform)

    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)
    test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=2)

    # Modele ResNet-18
    model = models.resnet18(weights="DEFAULT")
    for name, param in model.named_parameters():
        if "layer4" in name or "layer3" in name or "fc" in name:
            param.requires_grad = True
        else:
            param.requires_grad = False

    in_features = model.fc.in_features
    model.fc = nn.Sequential(
        nn.Linear(in_features, 256),
        nn.ReLU(),
        nn.Dropout(0.3),
        nn.Linear(256, 1),
    )
    model = model.to(device)

    trainable_params = [p for p in model.parameters() if p.requires_grad]
    optimizer = optim.Adam(trainable_params, lr=LR, weight_decay=WEIGHT_DECAY)
    criterion = nn.MSELoss()
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)

    best_loss = float("inf")
    patience_counter = 0
    history = []

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
        history.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss})

        if epoch % 5 == 0 or epoch == 1:
            print(f"  Fold {fold} epoch {epoch:3d} | train MSE = {train_loss:.4f} | val MSE = {val_loss:.4f}")

        if val_loss < best_loss:
            best_loss = val_loss
            torch.save(model.state_dict(), f"outputs/models/resnet18_fold{fold}.pt")
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= PATIENCE:
                print(f"  Fold {fold} early stopping at epoch {epoch}")
                break

    # Evaluation sur le test de ce fold
    model.load_state_dict(torch.load(f"outputs/models/resnet18_fold{fold}.pt"))
    model.eval()
    preds, trues = [], []
    with torch.no_grad():
        for x, y in test_loader:
            x = x.to(device)
            pred = model(x).squeeze(-1).cpu().numpy()
            preds.extend(pred)
            trues.extend(y.numpy())

    from sklearn.metrics import r2_score, mean_squared_error
    preds = np.array(preds); trues = np.array(trues)
    rmse = float(np.sqrt(mean_squared_error(trues, preds)))
    r2 = float(r2_score(trues, preds))
    print(f"  Fold {fold} -> R² = {r2:.3f}, RMSE = {rmse:.1f}")

    return preds, trues, test_df["cluster"].values, test_df["dept_id"].values, r2, rmse, history


def main():
    import torch
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import r2_score, mean_squared_error

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device : {device}")

    os.makedirs("outputs/tables", exist_ok=True)
    os.makedirs("outputs/models", exist_ok=True)
    os.makedirs("outputs/figures", exist_ok=True)

    # 1. Charger les labels depuis clusters_gps.csv
    print("=== 1. Chargement des donnees ===")
    if not os.path.exists(CLUSTERS_CSV):
        sys.exit(f"Fichier introuvable : {CLUSTERS_CSV}")
    labels = pd.read_csv(CLUSTERS_CSV)
    labels["cluster"] = labels["cluster"].astype(int)
    labels["dept_id"] = labels["dept_id"].astype(str)
    print(f"Clusters charges : {len(labels)}")

    # Verifier quels patches existent
    labels["patch_exists"] = labels["cluster"].apply(
        lambda c: load_patch(c, PATCH_DIR) is not None
    )
    labels = labels[labels["patch_exists"]].reset_index(drop=True)
    print(f"Clusters avec patch : {len(labels)}")

    if len(labels) < 10:
        sys.exit("Pas assez de patches")

    # 2. Group 5-fold par departement
    print(f"\n=== 2. Group {N_FOLDS}-fold par departement ===")
    gkf = GroupKFold(n_splits=N_FOLDS)
    groups = labels["dept_id"].values

    all_preds = np.empty(len(labels))
    all_trues = np.empty(len(labels))
    all_clusters = np.empty(len(labels))
    fold_r2 = []
    fold_rmse = []
    all_history = []

    for fold, (train_idx, test_idx) in enumerate(gkf.split(labels, groups=groups), 1):
        depts_test = labels.iloc[test_idx]["dept_id"].unique()
        print(f"\n--- Fold {fold}/{N_FOLDS} ({len(depts_test)} depts en test) ---")
        preds, trues, clusters, depts, r2, rmse, history = train_one_fold(
            train_idx, test_idx, labels, PATCH_DIR, device, fold
        )
        # Stocker les predictions OOF (out-of-fold)
        all_preds[test_idx] = preds
        all_trues[test_idx] = trues
        all_clusters[test_idx] = clusters
        fold_r2.append(r2)
        fold_rmse.append(rmse)
        all_history.append(history)

        # Nettoyer le checkpoint de ce fold
        ckpt = f"outputs/models/resnet18_fold{fold}.pt"
        if os.path.exists(ckpt):
            os.remove(ckpt)

    # 3. Metriques CV
    print("\n=== 3. Resultats de la validation croisee ===")
    r2_mean = float(np.mean(fold_r2))
    r2_std = float(np.std(fold_r2))
    rmse_mean = float(np.mean(fold_rmse))
    rmse_std = float(np.std(fold_rmse))

    print(f"R²  = {r2_mean:.3f} +/- {r2_std:.3f}")
    print(f"RMSE = {rmse_mean:.0f} +/- {rmse_std:.0f}")

    # R² global sur toutes les OOF predictions
    r2_global = float(r2_score(all_trues, all_preds))
    rmse_global = float(np.sqrt(mean_squared_error(all_trues, all_preds)))
    print(f"R² global (OOF) = {r2_global:.3f}, RMSE global (OOF) = {rmse_global:.0f}")

    # 4. Sauvegarder metriques
    pd.DataFrame([
        {"metric": "r2_mean_cv", "value": r2_mean},
        {"metric": "r2_std_cv", "value": r2_std},
        {"metric": "r2_global_oof", "value": r2_global},
        {"metric": "rmse_mean_cv", "value": rmse_mean},
        {"metric": "rmse_std_cv", "value": rmse_std},
        {"metric": "rmse_global_oof", "value": rmse_global},
        {"metric": "n_clusters", "value": len(labels)},
        {"metric": "n_folds", "value": N_FOLDS},
    ]).to_csv(OUT_METRICS, index=False)

    # 5. Figure : pred vs true (OOF)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(all_trues, all_preds, alpha=0.7, edgecolor="k")
    lo, hi = min(all_trues.min(), all_preds.min()), max(all_trues.max(), all_preds.max())
    ax.plot([lo, hi], [lo, hi], "r--", linewidth=1)
    ax.set_xlabel("True wealth (DHS)")
    ax.set_ylabel("Predicted wealth (ResNet-18)")
    ax.set_title(f"Group 5-fold CV: R² = {r2_global:.3f}, RMSE = {rmse_global:.0f}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG_PRED, dpi=200)
    plt.close(fig)

    # 6. Figure : boxplot des R² par fold
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.boxplot(fold_r2, vert=True, widths=0.4)
    ax.scatter([1] * len(fold_r2), fold_r2, color="red", zorder=5, s=40)
    ax.set_ylabel("R²")
    ax.set_title(f"R² per fold (mean = {r2_mean:.3f}, sd = {r2_std:.3f})")
    ax.set_xticks([1])
    ax.set_xticklabels([f"{N_FOLDS}-fold CV"])
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG_CV, dpi=200)
    plt.close(fig)

    # 7. Courbes d'apprentissage (moyennees sur tous les folds)
    max_ep = max(len(h) for h in all_history)
    epochs = np.arange(1, max_ep + 1)
    train_losses = np.full((len(all_history), max_ep), np.nan)
    val_losses = np.full((len(all_history), max_ep), np.nan)
    for i, h in enumerate(all_history):
        for entry in h:
            train_losses[i, entry["epoch"] - 1] = entry["train_loss"]
            val_losses[i, entry["epoch"] - 1] = entry["val_loss"]

    train_mean = np.nanmean(train_losses, axis=0)
    train_std = np.nanstd(train_losses, axis=0)
    val_mean = np.nanmean(val_losses, axis=0)
    val_std = np.nanstd(val_losses, axis=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(epochs, train_mean, label="Train (mean)", color="blue")
    ax.fill_between(epochs, train_mean - train_std, train_mean + train_std, alpha=0.2, color="blue")
    ax.plot(epochs, val_mean, label="Validation (mean)", color="orange")
    ax.fill_between(epochs, val_mean - val_std, val_mean + val_std, alpha=0.2, color="orange")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("MSE")
    ax.set_title("Learning curves (mean across 5 folds)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG_CURVES, dpi=200)
    plt.close(fig)

    # Sauvegarder l'historique d'apprentissage
    history_rows = []
    for fold_idx, h in enumerate(all_history):
        for entry in h:
            history_rows.append({
                "fold": fold_idx + 1,
                "epoch": entry["epoch"],
                "train_loss": entry["train_loss"],
                "val_loss": entry["val_loss"],
            })
    pd.DataFrame(history_rows).to_csv(OUT_HISTORY, index=False)

    # 8. Sauvegarder predictions OOF pour le FH
    preds_df = pd.DataFrame({
        "cluster": all_clusters.astype(int),
        "wealth_true": all_trues,
        "wealth_pred": all_preds,
    })
    preds_df.to_csv(OUT_PREDS, index=False)

    # 9. Predire sur TOUS les clusters (pour le FH)
    print("\n=== 9. Entrainement final sur toutes les donnees ===")
    # Reinitialiser le modele et entrainer sur tout
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
    import torchvision.models as models
    import torchvision.transforms as T

    class LandsatDatasetFull(Dataset):
        def __init__(self, df, patch_dir, transform=None):
            self.df = df.reset_index(drop=True)
            self.patch_dir = patch_dir
            self.transform = transform

        def __len__(self):
            return len(self.df)

        def __getitem__(self, idx):
            row = self.df.iloc[idx]
            patch = load_patch(row["cluster"], self.patch_dir)
            patch = torch.from_numpy(patch).permute(2, 0, 1).float()
            if self.transform:
                patch = self.transform(patch)
            return patch, torch.tensor(row["wealth_mean"], dtype=torch.float32)

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

    full_ds = LandsatDatasetFull(labels, PATCH_DIR, transform=train_transform)
    full_loader = DataLoader(full_ds, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)

    model = models.resnet18(weights="DEFAULT")
    for name, param in model.named_parameters():
        if "layer4" in name or "layer3" in name or "fc" in name:
            param.requires_grad = True
        else:
            param.requires_grad = False

    in_features = model.fc.in_features
    model.fc = nn.Sequential(
        nn.Linear(in_features, 256),
        nn.ReLU(),
        nn.Dropout(0.3),
        nn.Linear(256, 1),
    )
    model = model.to(device)

    trainable_params = [p for p in model.parameters() if p.requires_grad]
    optimizer = optim.Adam(trainable_params, lr=LR, weight_decay=WEIGHT_DECAY)
    criterion = nn.MSELoss()

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        train_loss = 0.0
        for x, y in full_loader:
            x, y = x.to(device), y.to(device)
            optimizer.zero_grad()
            pred = model(x).squeeze(-1)
            loss = criterion(pred, y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * x.size(0)
        train_loss /= len(full_ds)
        if epoch % 10 == 0 or epoch == 1:
            print(f"  Full epoch {epoch:3d} | train MSE = {train_loss:.4f}")

    torch.save(model.state_dict(), OUT_MODEL)
    print(f"Modele final : {OUT_MODEL}")

    # Predire sur tous les clusters
    model.eval()
    eval_ds = LandsatDatasetFull(labels, PATCH_DIR, transform=test_transform)
    eval_loader = DataLoader(eval_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=2)
    all_full_preds = []
    with torch.no_grad():
        for x, _ in eval_loader:
            x = x.to(device)
            pred = model(x).squeeze(-1).cpu().numpy()
            all_full_preds.extend(pred)

    pd.DataFrame({
        "cluster": labels["cluster"].values,
        "cnn_prediction": all_full_preds,
    }).to_csv("data/processed/cnn_predictions_cluster.csv", index=False)
    print("Predictions tous clusters : data/processed/cnn_predictions_cluster.csv")

    print("\nTermine.")
    print(f"Metriques CV : {OUT_METRICS}")
    print(f"Predictions OOF : {OUT_PREDS}")
    print(f"Figures : {OUT_FIG_PRED}, {OUT_FIG_CURVES}, {OUT_FIG_CV}")
    print(f"Historique : {OUT_HISTORY}")


if __name__ == "__main__":
    main()
