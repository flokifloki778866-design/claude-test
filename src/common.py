"""Shared loaders and plotting defaults for the CPEC satellite files.

All five GEE exports carry two junk columns: `system:index` (an Earth Engine row id)
and `.geo` (always an empty MultiPoint, i.e. no coordinates were exported).
Both are dropped on load.
"""
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CROSSWALK = ROOT / "data" / "crosswalk"
OUT = ROOT / "output"

FILES = {
    "ntl_district": "1_ntl_district_annual.csv",
    "ntl_projects": "2_ntl_rings_projects.csv",
    "ntl_roads": "3_ntl_rings_roads.csv",
    "builtup_district": "4_builtup_district.csv",
    "no2_monthly": "6_no2_monthly.csv",
}

RING_LABELS = {0: "0–10 km", 10: "10–25 km", 25: "25–50 km", 50: "50–100 km"}

# Categorical slots (fixed order) plus a recessive gray for context series.
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
GRAY = "#a3a29c"
INK, INK2 = "#0b0b0b", "#52514e"


def load(key: str) -> pd.DataFrame:
    d = pd.read_csv(RAW / FILES[key])
    return d.drop(columns=["system:index", ".geo"], errors="ignore")


def style():
    plt.rcParams.update({
        "figure.dpi": 130,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": GRAY,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.grid": True,
        "grid.color": "#e6e5e0",
        "grid.linewidth": 0.6,
        "lines.linewidth": 2,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "legend.frameon": False,
    })
