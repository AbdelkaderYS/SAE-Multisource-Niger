# ===========================================================================
# gradcam.py
#
# Visualise les zones du patch Sentinel-2 qui activent le plus le neurone
# de sortie du CNN entraine. La methode Grad-CAM (Gradient-weighted Class
# Activation Mapping, Selvaraju et al. 2017) est appliquee sur la derniere
# couche convolutive.
#
# Sorties : outputs/figures/04_gradcam_<statut>_<idx>.png
#   ou statut = "poor" / "rich" (2 clusters les plus pauvres et 2 plus riches)
# ===========================================================================

import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cv2
except ImportError:
    sys.exit("OpenCV requis : pip install opencv-python")


PATCH_DIR = "python/data/patches"
TABLE_CSV = "python/data/clusters_wealth_gps.csv"
MODEL_PT  = "python/data/cnn_model.pt"
OUT_DIR   = "outputs/figures"

# Reprendre la meme architecture que dans entrainer_cnn.py
from entrainer_cnn import SmallCNN, ClusterPatchDataset  # noqa: E402


class GradCAM:
    def __init__(self, model, target_layer):
        self.model = model
        self.target_layer = target_layer
        self.gradients = None
        self.activations = None
        target_layer.register_forward_hook(self._save_activation)
        target_layer.register_full_backward_hook(self._save_gradient)

    def _save_activation(self, module, inp, out):
        self.activations = out.detach()

    def _save_gradient(self, module, grad_in, grad_out):
        self.gradients = grad_out[0].detach()

    def __call__(self, x):
        self.model.zero_grad()
        out = self.model(x)
        out.sum().backward()
        weights = self.gradients.mean(dim=(2, 3), keepdim=True)
        cam = (weights * self.activations).sum(dim=1, keepdim=True)
        cam = F.relu(cam)
        cam = F.interpolate(cam, size=x.shape[2:], mode="bilinear", align_corners=False)
        cam = cam.squeeze().cpu().numpy()
        cam = (cam - cam.min()) / (cam.max() - cam.min() + 1e-8)
        return cam


def overlay(rgb, cam):
    heatmap = plt.cm.jet(cam)[:, :, :3]
    overlay = 0.5 * rgb + 0.5 * heatmap
    overlay = np.clip(overlay, 0, 1)
    return overlay


def main():
    if not os.path.exists(MODEL_PT):
        sys.exit(f"Modele introuvable : {MODEL_PT}\nLancer d'abord entrainer_cnn.py")
    if not os.path.exists(TABLE_CSV):
        sys.exit(f"Table introuvable : {TABLE_CSV}\nLancer d'abord preparer_clusters.py")
    if not os.path.isdir(PATCH_DIR) or len(os.listdir(PATCH_DIR)) == 0:
        sys.exit(f"Patches introuvables dans {PATCH_DIR}.")

    df = pd.read_csv(TABLE_CSV)
    df["cluster"] = df["cluster"].astype(int)

    # Charger le modele
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = SmallCNN().to(device)
    model.load_state_dict(torch.load(MODEL_PT, map_location=device))
    model.eval()

    # Couche cible : derniere conv du features
    target_layer = model.features[6]
    gradcam = GradCAM(model, target_layer)

    # Selectionner 2 clusters les plus pauvres et 2 plus riches
    df_sorted = df.sort_values("wealth_mean")
    poor_clusters = df_sorted.head(2)["cluster"].tolist()
    rich_clusters = df_sorted.tail(2)["cluster"].tolist()

    # Utiliser les memes stats que lors de l'entrainement
    train_df, _ = None, None
    from sklearn.model_selection import train_test_split
    train_df, _ = train_test_split(df, test_size=0.20, random_state=42)
    ds_stats = ClusterPatchDataset(train_df, PATCH_DIR)
    mean, std = ds_stats.mean, ds_stats.std

    os.makedirs(OUT_DIR, exist_ok=True)

    def visualize(clusters, statut):
        for idx, cluster in enumerate(clusters, start=1):
            arr = np.load(os.path.join(PATCH_DIR, f"{int(cluster)}.npy"))
            # RGB composite : B4 (rouge), B3 (vert), B2 (bleu)
            rgb = np.stack([arr[2], arr[3], arr[1]], axis=-1)
            rgb = (rgb - rgb.min()) / (rgb.max() - rgb.min() + 1e-8)
            x = (arr - mean[:, None, None]) / std[:, None, None]
            x_t = torch.from_numpy(x).float().unsqueeze(0).to(device)
            cam = gradcam(x_t)
            over = overlay(rgb, cam)
            fig, ax = plt.subplots(1, 2, figsize=(10, 5))
            ax[0].imshow(rgb); ax[0].set_title(f"Patch Sentinel-2\nCluster {cluster}")
            ax[1].imshow(over); ax[1].set_title("Grad-CAM")
            for a in ax: a.axis("off")
            fig.tight_layout()
            out = os.path.join(OUT_DIR, f"04_gradcam_{statut}_{idx}.png")
            fig.savefig(out, dpi=200)
            plt.close(fig)
            print(f"Grad-CAM sauvegarde : {out}")

    visualize(poor_clusters, "poor")
    visualize(rich_clusters, "rich")
    print("Termine.")


if __name__ == "__main__":
    main()
