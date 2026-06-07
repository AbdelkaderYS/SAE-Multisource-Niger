# ===========================================================================
# entrainer_cnn.py
#
# Entraine un CNN simple a 3 couches de convolution pour predire la richesse
# moyenne d'un cluster DHS a partir d'un patch Sentinel-2 (4 bandes, 64x64).
#
# Architecture (volontairement simple, ~5k parametres) :
#   Input  (B, 4, 64, 64)
#     -> Conv2D(32, 3) + ReLU + MaxPool(2)   -> (B, 32, 31, 31)
#     -> Conv2D(64, 3) + ReLU + MaxPool(2)   -> (B, 64, 14, 14)
#     -> Conv2D(128, 3) + ReLU + MaxPool(2)  -> (B, 128, 6, 6)
#     -> Flatten + Linear(256) + ReLU + Dropout(0.5) + Linear(1)
#
# Sorties :
#   outputs/tables/cnn_metrics.csv
#   outputs/tables/cnn_predictions.csv
#   outputs/figures/03_cnn_pred_vs_true.png
#   python/data/cnn_model.pt
# ===========================================================================

import os
import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import torch
    import torch.nn as nn
    from torch.utils.data import Dataset, DataLoader
except ImportError:
    sys.exit("PyTorch requis : pip install torch torchvision")


PATCH_DIR = "python/data/patches"
TABLE_CSV = "python/data/clusters_wealth_gps.csv"
OUT_METRICS = "outputs/tables/cnn_metrics.csv"
OUT_PREDS   = "outputs/tables/cnn_predictions.csv"
OUT_FIG     = "outputs/figures/03_cnn_pred_vs_true.png"
OUT_MODEL   = "python/data/cnn_model.pt"


class ClusterPatchDataset(Dataset):
    def __init__(self, df, patch_dir, stats=None):
        self.df = df.reset_index(drop=True)
        self.patch_dir = patch_dir
        # Normalisation par canal (mean/std sur le train)
        if stats is None:
            self.mean, self.std = self._compute_stats()
        else:
            self.mean, self.std = stats

    def _compute_stats(self):
        all_pixels = []
        for cluster in self.df["cluster"]:
            arr = np.load(os.path.join(self.patch_dir, f"{int(cluster)}.npy"))
            all_pixels.append(arr.reshape(4, -1).T)
        all_pixels = np.concatenate(all_pixels, axis=0)
        return all_pixels.mean(axis=0), all_pixels.std(axis=0) + 1e-8

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        arr = np.load(os.path.join(self.patch_dir, f"{int(row['cluster'])}.npy"))
        # arr : (4, 64, 64) float
        arr = (arr - self.mean[:, None, None]) / self.std[:, None, None]
        return torch.from_numpy(arr).float(), torch.tensor(row["wealth_mean"], dtype=torch.float32)


class SmallCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(4,  32, kernel_size=3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
            nn.Conv2d(32, 64, kernel_size=3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
            nn.Conv2d(64, 128, kernel_size=3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
        )
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.Linear(128 * 8 * 8, 256), nn.ReLU(), nn.Dropout(0.5),
            nn.Linear(256, 1),
        )

    def forward(self, x):
        return self.head(self.features(x)).squeeze(-1)


def train_one_epoch(model, loader, opt, loss_fn, device):
    model.train()
    total = 0.0
    for x, y in loader:
        x, y = x.to(device), y.to(device)
        opt.zero_grad()
        pred = model(x)
        loss = loss_fn(pred, y)
        loss.backward()
        opt.step()
        total += loss.item() * x.size(0)
    return total / len(loader.dataset)


def evaluate(model, loader, device):
    model.eval()
    preds, trues = [], []
    with torch.no_grad():
        for x, y in loader:
            x = x.to(device)
            p = model(x).cpu().numpy()
            preds.append(p)
            trues.append(y.numpy())
    return np.concatenate(preds), np.concatenate(trues)


def main():
    if not os.path.exists(TABLE_CSV):
        sys.exit(f"Table introuvable : {TABLE_CSV}\nLancer d'abord preparer_clusters.py")
    if not os.path.isdir(PATCH_DIR) or len(os.listdir(PATCH_DIR)) == 0:
        sys.exit(f"Aucun patch Sentinel-2 dans {PATCH_DIR}.\n"
                 "Televerser via Google Earth Engine (voir preparer_clusters.py).")

    df = pd.read_csv(TABLE_CSV)
    df = df[df["cluster"].notna()].copy()
    df["cluster"] = df["cluster"].astype(int)
    print(f"Clusters disponibles : {len(df)}")

    # Split 80/20
    train_df, test_df = train_test_split(df, test_size=0.20, random_state=42)
    print(f"Train : {len(train_df)} | Test : {len(test_df)}")

    train_ds = ClusterPatchDataset(train_df, PATCH_DIR)
    test_ds  = ClusterPatchDataset(test_df,  PATCH_DIR, stats=(train_ds.mean, train_ds.std))

    train_loader = DataLoader(train_ds, batch_size=32, shuffle=True)
    test_loader  = DataLoader(test_ds,  batch_size=32, shuffle=False)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device : {device}")

    model = SmallCNN().to(device)
    opt   = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
    loss_fn = nn.MSELoss()

    n_epochs = 20
    for epoch in range(1, n_epochs + 1):
        loss = train_one_epoch(model, train_loader, opt, loss_fn, device)
        if epoch % 5 == 0 or epoch == 1:
            print(f"Epoch {epoch:02d} | train MSE = {loss:.4f}")

    # Evaluation finale
    test_preds, test_trues = evaluate(model, test_loader, device)
    rmse = float(np.sqrt(mean_squared_error(test_trues, test_preds)))
    r2   = float(r2_score(test_trues, test_preds))
    print(f"\nTest RMSE : {rmse:.3f}")
    print(f"Test R2   : {r2:.3f}")

    # Sauvegardes
    os.makedirs("outputs/tables", exist_ok=True)
    os.makedirs("outputs/figures", exist_ok=True)
    os.makedirs("python/data", exist_ok=True)

    pd.DataFrame([{"metric": "rmse_test", "value": rmse},
                  {"metric": "r2_test",   "value": r2}]).to_csv(OUT_METRICS, index=False)
    pd.DataFrame({"cluster": test_df["cluster"].values,
                  "wealth_true": test_trues,
                  "wealth_pred": test_preds}).to_csv(OUT_PREDS, index=False)
    torch.save(model.state_dict(), OUT_MODEL)

    # Figure
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(test_trues, test_preds, alpha=0.7, edgecolor="k")
    lo = min(test_trues.min(), test_preds.min())
    hi = max(test_trues.max(), test_preds.max())
    ax.plot([lo, hi], [lo, hi], "r--", linewidth=1)
    ax.set_xlabel("Richesse vraie (DHS)")
    ax.set_ylabel("Richesse predite (CNN)")
    ax.set_title(f"CNN sur Sentinel-2 : Test RMSE = {rmse:.3f}, R2 = {r2:.3f}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=200)
    print(f"Figure sauvegardee : {OUT_FIG}")
    print(f"Metriques : {OUT_METRICS}")
    print(f"Modele    : {OUT_MODEL}")


if __name__ == "__main__":
    main()
