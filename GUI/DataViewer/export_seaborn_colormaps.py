# export_seaborn_colormaps.py
"""
Export selected seaborn/matplotlib colormaps for MATLAB App Designer.

Outputs:
    resources/colormaps/seaborn_color_palettes.mat
    resources/colormaps/icons/<colormap_name>.png

All colormaps are saved in MATLAB-standard Nx3 format (here: 64x3).
"""

from pathlib import Path
import re

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.io import savemat


# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------
N_COLORS = 64

# Change this if you want another output location
OUTPUT_DIR = Path(__file__).resolve().parent / "colormaps"
ICON_DIR = OUTPUT_DIR / "icons"
MAT_FILE = OUTPUT_DIR / "seaborn_color_palettes.mat"


# -------------------------------------------------------------------------
# Palette lists
# -------------------------------------------------------------------------
PERCEPTUALLY_UNIFORM = [
    "rocket",
    "mako",
    "flare",
    "crest",
    "viridis",
    "magma",
    "inferno",
    "plasma",
    "cividis",
]

DIVERGENT = [
    "vlag",
    "icefire",
    "coolwarm",
    "Spectral",
    "RdBu_r",
    "RdYlBu",
    "PiYG",
    "seismic",
]

PALETTES = PERCEPTUALLY_UNIFORM + DIVERGENT


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
def make_matlab_safe_name(name: str) -> str:
    """Return a valid MATLAB variable/struct field name."""
    safe = re.sub(r"\W+", "_", name)
    if not re.match(r"^[A-Za-z]", safe):
        safe = f"cmap_{safe}"
    return safe


def sample_colormap(name: str, n_colors: int) -> np.ndarray:
    """
    Sample a seaborn or matplotlib colormap as an RGB array.

    Returns:
        n_colors x 3 float64 array in [0, 1]
    """
    try:
        cmap = sns.color_palette(name, as_cmap=True)
        rgb = np.asarray(cmap(np.linspace(0, 1, n_colors)))[:, :3]
    except Exception:
        cmap = mpl.colormaps.get_cmap(name)
        rgb = np.asarray(cmap(np.linspace(0, 1, n_colors)))[:, :3]

    rgb = np.asarray(rgb, dtype=np.float64)
    rgb = np.clip(rgb, 0.0, 1.0)
    return rgb


def save_colormap_icon(cmap_rgb: np.ndarray, icon_path: Path) -> None:
    """Save a horizontal gradient PNG preview for one colormap."""
    icon_path.parent.mkdir(parents=True, exist_ok=True)
    height = 18
    gradient = np.tile(cmap_rgb[np.newaxis, :, :], (height, 1, 1))
    plt.imsave(icon_path, gradient)


# -------------------------------------------------------------------------
# Main export
# -------------------------------------------------------------------------
def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ICON_DIR.mkdir(parents=True, exist_ok=True)

    mat_payload = {}
    cmap_names = []
    cmap_display_names = []
    cmap_icon_files = []
    cmap_groups = []

    for name in PALETTES:
        safe_name = make_matlab_safe_name(name)
        cmap_rgb = sample_colormap(name, N_COLORS)   # Always Nx3

        icon_path = ICON_DIR / f"{safe_name}.png"
        save_colormap_icon(cmap_rgb, icon_path)

        mat_payload[safe_name] = cmap_rgb
        cmap_names.append(safe_name)
        cmap_display_names.append(name)
        cmap_icon_files.append(str(icon_path))
        cmap_groups.append(
            "perceptually_uniform"
            if name in PERCEPTUALLY_UNIFORM
            else "divergent"
        )

        print(f"Exported {name:12s} -> {safe_name:12s}")

    # Metadata for MATLAB-side setup
    mat_payload["cmapNames"] = np.asarray(cmap_names, dtype=object)
    mat_payload["cmapDisplayNames"] = np.asarray(cmap_display_names, dtype=object)
    mat_payload["cmapIconFiles"] = np.asarray(cmap_icon_files, dtype=object)
    mat_payload["cmapGroups"] = np.asarray(cmap_groups, dtype=object)
    mat_payload["nColors"] = np.asarray([[N_COLORS]], dtype=np.float64)

    savemat(MAT_FILE, mat_payload, do_compression=True)

    print(f"\nSaved MAT file: {MAT_FILE}")
    print(f"Saved icons to: {ICON_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()