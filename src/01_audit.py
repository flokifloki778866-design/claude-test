"""Week-1 data audit: dictionary, coverage, duplicates, sensor break, raw trend plots,
and the project-UID crosswalk template.

Run from the repo root:  python src/01_audit.py
Outputs go to output/audit/ and data/crosswalk/.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import (CROSSWALK, OUT, RAW, FILES, RING_LABELS, C1, C2, GRAY, INK2,
                    load, style)

AUD = OUT / "audit"
AUD.mkdir(parents=True, exist_ok=True)
CROSSWALK.mkdir(parents=True, exist_ok=True)
style()

# Ring bands as a sequential ramp (near = dark) - ring distance is ordered, not categorical.
RING_COLORS = {0: "#0d3d7a", 10: "#2a78d6", 25: "#7fb0ea", 50: "#c4dbf6"}

KEYS = {
    "ntl_district": ["district", "province", "year"],
    "ntl_projects": ["uid", "ring_lo_km", "year"],
    "ntl_roads": ["uid", "ring_lo_km", "year"],
    "builtup_district": ["district", "province", "year"],
    "no2_monthly": ["name", "ring", "date"],
}


# 1. Data dictionary and coverage -------------------------------------------------------
dict_rows, cov_rows = [], []
data = {k: load(k) for k in FILES}
for k, d in data.items():
    raw = pd.read_csv(RAW / FILES[k])
    geo_empty = (raw[".geo"] == '{"type":"MultiPoint","coordinates":[]}').all()
    for c in d.columns:
        s = d[c]
        num = pd.api.types.is_numeric_dtype(s)
        dict_rows.append({
            "file": FILES[k], "column": c, "dtype": str(s.dtype),
            "n_unique": s.nunique(), "n_missing": int(s.isna().sum()),
            "min": s.min() if num else "", "max": s.max() if num else "",
            "examples": "; ".join(map(str, s.drop_duplicates().head(4).tolist())),
        })
    unit = KEYS[k][:-1]
    cov_rows.append({
        "file": FILES[k], "rows": len(d),
        "units": d.drop_duplicates(unit).shape[0],
        "periods": d[KEYS[k][-1]].nunique(),
        "first": d[KEYS[k][-1]].min(), "last": d[KEYS[k][-1]].max(),
        "balanced": d.groupby(unit).size().nunique() == 1,
        "duplicate_keys": int(d.duplicated(KEYS[k]).sum()),
        "missing_cells": int(d.isna().sum().sum()),
        "coordinates_exported": not geo_empty,
    })
pd.DataFrame(dict_rows).to_csv(AUD / "data_dictionary.csv", index=False)
cov = pd.DataFrame(cov_rows)
cov.to_csv(AUD / "coverage.csv", index=False)
print(cov.to_string(index=False))


# 2. DMSP -> VIIRS break ------------------------------------------------------------------
def cell_id(d):
    if "ring_lo_km" in d:
        return d["uid"].astype(str) + "_" + d["ring_lo_km"].astype(str)
    return d["district"] + "_" + d["province"]


brk = []
for k in ["ntl_district", "ntl_projects", "ntl_roads"]:
    d = data[k].assign(cell=lambda x: cell_id(x))
    w = d.pivot_table(index="cell", columns="year", values="ntl_mean")
    lw = np.log(w + 0.01)
    row = {"file": FILES[k], "cells": len(w),
           "dmsp_years": f"{d[d.sensor == 'DMSP'].year.min()}-{d[d.sensor == 'DMSP'].year.max()}",
           "viirs_years": f"{d[d.sensor == 'VIIRS'].year.min()}-{d[d.sensor == 'VIIRS'].year.max()}",
           "overlap_years": len(set(d[d.sensor == "DMSP"].year) & set(d[d.sensor == "VIIRS"].year)),
           "dmsp_share_zero": round((d[d.sensor == "DMSP"].ntl_mean == 0).mean(), 3),
           "dmsp_max": round(d[d.sensor == "DMSP"].ntl_mean.max(), 2),
           "dmsp_cells_at_63": int((d[d.sensor == "DMSP"].ntl_mean >= 62.5).sum())}
    for a, b in [(2011, 2012), (2012, 2013), (2013, 2014), (2014, 2015), (2015, 2016)]:
        row[f"median_dlog_{a}_{b}"] = round(float(np.nanmedian(lw[b] - lw[a])), 3)
    row["corr_log_2013_2014"] = round(float(np.corrcoef(lw[2013], lw[2014])[0, 1]), 3)
    row["corr_log_2012_2013"] = round(float(np.corrcoef(lw[2012], lw[2013])[0, 1]), 3)
    brk.append(row)
brk = pd.DataFrame(brk)
brk.to_csv(AUD / "sensor_break.csv", index=False)
print("\n", brk.T.to_string())

fig, ax = plt.subplots(figsize=(7, 3.6))
for k, col, lab in [("ntl_district", C1, "Districts (119)"),
                    ("ntl_projects", C2, "Project rings (43 × 4)"),
                    ("ntl_roads", GRAY, "Road rings (6 × 4)")]:
    m = data[k].groupby("year").ntl_mean.median()
    ax.plot(m.index, m.values, color=col, label=lab, marker="o", markersize=3)
ax.axvline(2013.5, color=INK2, lw=1, ls="--")
ax.text(2013.7, ax.get_ylim()[1] * 0.9, "DMSP → VIIRS\n(no overlap year)", color=INK2, fontsize=8, va="top")
ax.set_yscale("log")
ax.set_ylabel("Median ntl_mean (log scale)")
ax.set_title("Nightlights are raw DMSP DN, then raw VIIRS radiance: a ~9× level break in 2014")
ax.legend(loc="lower left")
fig.tight_layout()
fig.savefig(AUD / "fig1_sensor_break.png")
plt.close(fig)


# 3. Raw trends by ring -------------------------------------------------------------------
def ring_panel(ax, d, title):
    for lo, lab in RING_LABELS.items():
        m = d[d.ring_lo_km == lo].groupby("year").ntl_mean.median()
        ax.plot(m.index, m.values, color=RING_COLORS[lo], label=lab, lw=1.6)
    ax.axvline(2013.5, color=INK2, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_title(title)


fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
ring_panel(axes[0], data["ntl_projects"], "Project rings: median across 43 UIDs")
ring_panel(axes[1], data["ntl_roads"], "Road rings: median across 6 CPEC roads")
axes[0].legend(title="Ring", fontsize=8)
fig.tight_layout()
fig.savefig(AUD / "fig2_rings_raw.png")
plt.close(fig)

roads = data["ntl_roads"]
fig, axes = plt.subplots(2, 3, figsize=(10, 5.4), sharex=True)
for ax, uid in zip(axes.flat, sorted(roads.uid.unique())):
    ring_panel(ax, roads[roads.uid == uid], uid)
axes.flat[0].legend(title="Ring", fontsize=7)
fig.suptitle("Road rings by route (no placebo roads M-1/M-2/M-3/N-5 in the export)", fontweight="bold")
fig.tight_layout()
fig.savefig(AUD / "fig3_roads_by_route.png")
plt.close(fig)

proj = data["ntl_projects"]
uids = sorted(proj.uid.unique())
fig, axes = plt.subplots(6, 8, figsize=(16, 11), sharex=True)
for ax in axes.flat:
    ax.set_visible(False)
for ax, uid in zip(axes.flat, uids):
    ax.set_visible(True)
    ring_panel(ax, proj[proj.uid == uid], str(uid))
    ax.tick_params(labelsize=6)
    ax.title.set_fontsize(8)
axes.flat[0].legend(fontsize=6)
fig.suptitle("Project rings by UID, raw ntl_mean (log). Labels and treatment dates still unknown.",
             fontweight="bold")
fig.tight_layout()
fig.savefig(AUD / "fig4_projects_by_uid.png")
plt.close(fig)


# 3b. Duplicate-location clusters -----------------------------------------------------------
# Many UIDs have identical ring series, i.e. they were buffered around the same point.
# Cluster UIDs whose full 4-ring x 33-year log series differ by < 0.02 on average.
lw = np.log(proj.pivot_table(index="uid", columns=["ring_lo_km", "year"], values="ntl_mean") + 0.01)
dist = np.abs(lw.values[:, None, :] - lw.values[None, :, :]).mean(axis=2)
site = {}
for i, uid in enumerate(lw.index):
    if uid in site:
        continue
    sid = f"S{len(set(site.values())) + 1:02d}"
    for j in np.where(dist[i] < 0.02)[0]:
        site.setdefault(lw.index[j], sid)
sites = pd.Series(site, name="site_id").rename_axis("uid").reset_index()
sites["n_uids_at_site"] = sites.groupby("site_id").uid.transform("size")
sites["exact_duplicate_series"] = sites.groupby("site_id").uid.transform(
    lambda s: lw.loc[s].round(6).drop_duplicates().shape[0] == 1)
sites.to_csv(AUD / "project_site_clusters.csv", index=False)
print(f"\n{proj.uid.nunique()} UIDs -> {sites.site_id.nunique()} distinct locations; "
      f"largest cluster {sites.n_uids_at_site.max()} UIDs")


# 4. Built-up -----------------------------------------------------------------------------
b = data["builtup_district"]
bw = b.pivot_table(index=["district", "province"], columns="year", values="built_m2")
bu = pd.DataFrame({
    "national_km2": (bw.sum() / 1e6).round(0),
    "districts_declining_vs_prev_epoch": (bw.diff(axis=1) < 0).sum(),
    "median_5y_growth_pct": (bw.pct_change(axis=1) * 100).median().round(1),
})
bu.to_csv(AUD / "builtup_epochs.csv")
print("\n", bu.to_string())

fig, ax = plt.subplots(figsize=(6, 3.2))
ax.plot(bu.index, bu.national_km2, color=C1, marker="o")
ax.plot([2020, 2025], bu.national_km2.loc[[2020, 2025]], color="white", lw=2.2)
ax.plot([2020, 2025], bu.national_km2.loc[[2020, 2025]], color=C1, lw=2, ls=":")
ax.annotate("2025 epoch: check whether\nobserved or GHSL projection", (2025, bu.national_km2[2025]),
            xytext=(2008, bu.national_km2[2025] * 0.95), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=GRAY))
ax.set_ylabel("Built-up area, km²")
ax.set_title("Built-up area, sum of 119 districts (5-year epochs)")
fig.tight_layout()
fig.savefig(AUD / "fig5_builtup_national.png")
plt.close(fig)


# 5. District list for the LFS crosswalk ---------------------------------------------------
(data["ntl_district"][["district", "province"]].drop_duplicates()
 .sort_values(["province", "district"])
 .assign(lfs_unit_86="", note="")
 .to_csv(CROSSWALK / "district119_to_lfs86_TEMPLATE.csv", index=False))


# 6. Project UID crosswalk template -------------------------------------------------------
v = proj[proj.sensor == "VIIRS"]
summ = (v.pivot_table(index="uid", columns="ring_lo_km", values="ntl_mean", aggfunc="mean")
        .rename(columns=lambda c: f"viirs_mean_{RING_LABELS[c].replace(' ', '').replace('–', '_')}"))
tmpl = pd.DataFrame(index=pd.Index(uids, name="uid"))
for c in ["project_name", "aiddata_title", "cpec_tagged", "sector", "project_type",
          "capacity_mw_or_km", "lat", "lon", "announcement_date", "construction_start",
          "commissioning_date", "treatment_date_source", "treated_or_control", "notes"]:
    tmpl[c] = ""
tmpl = sites.set_index("uid")[["site_id", "n_uids_at_site"]].join(tmpl, how="right")
tmpl = tmpl.join(summ.round(3))
tmpl.to_csv(CROSSWALK / "project_uid_crosswalk_TEMPLATE.csv")
print(f"\nWrote crosswalk template for {len(uids)} UIDs")
