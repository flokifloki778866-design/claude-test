"""Week-2 NO2 diagnostics: core-minus-control series, pre-period inventory, seasonality,
and per-site deseasonalised pre/post contrasts with a fake-date placebo.

These are feasibility diagnostics, not causal estimates: 2-4 usable treated sites,
no meteorology, and ring radii/coordinates not in the export.

Run from the repo root:  python src/02_no2_diagnostics.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

from common import OUT, C1, C2, GRAY, INK2, load, style

NO2 = OUT / "no2"
NO2.mkdir(parents=True, exist_ok=True)
style()

n = load("no2_monthly")
n["date"] = pd.to_datetime(n["date"])
n["commissioned"] = pd.to_datetime(n["commissioned"])
n["umol"] = n["value"] * 1e6  # mol/m^2 -> umol/m^2

w = n.pivot_table(index=["name", "fuel", "mw", "commissioned", "date", "month", "year"],
                  columns="ring", values="umol").reset_index()
w["diff"] = w["core"] - w["control"]
w["logratio"] = np.log(w["core"] / w["control"])
w["cod_month"] = w["commissioned"].dt.to_period("M").dt.to_timestamp()
w["rel_month"] = ((w["date"].dt.year - w["cod_month"].dt.year) * 12
                  + (w["date"].dt.month - w["cod_month"].dt.month))
w["post"] = (w["rel_month"] >= 0).astype(int)
w.to_csv(NO2 / "no2_core_control_panel.csv", index=False)

SITES = ["Sahiwal", "PortQasim", "Hub", "EngroThar", "SSRL_Thar", "ThalNova",
         "Karot_hydro", "SukiKinari_hydro"]

# 1. Inventory and cross-site dependence -----------------------------------------------------
core = w.pivot_table(index="date", columns="name", values="core")
corr_core = core.corr()
corr_core.round(2).to_csv(NO2 / "core_ring_correlation.csv")

inv = (w.groupby(["name", "fuel", "mw", "commissioned"])
       .agg(pre_months=("post", lambda s: int((s == 0).sum())),
            post_months=("post", "sum"),
            core_mean=("core", "mean"), control_mean=("control", "mean"),
            diff_mean=("diff", "mean"))
       .round(2).reset_index().set_index("name").loc[SITES])
inv["max_core_corr_other_site"] = [corr_core.loc[s].drop(s).max().round(2) for s in inv.index]
inv["most_correlated_site"] = [corr_core.loc[s].drop(s).idxmax() for s in inv.index]


def usable(r):
    if r.pre_months == 0:
        return "no pre-period (COD before S5P series starts)"
    if r.max_core_corr_other_site > 0.95:
        return f"same footprint as {r.most_correlated_site}: not an independent site"
    if r.pre_months < 24:
        return "usable, but <2 pre-years (one seasonal cycle)"
    return "usable"


inv["assessment"] = inv.apply(usable, axis=1)
inv.to_csv(NO2 / "site_inventory.csv")
print(inv.to_string())


# 2. Per-site deseasonalised pre/post contrast, with fake-date placebo ----------------------
def contrast(d, y, post_col):
    m = smf.ols(f"{y} ~ {post_col} + C(month)", data=d).fit(cov_type="HAC", cov_kwds={"maxlags": 6})
    return m.params[post_col], m.bse[post_col]


rows = []
for s in SITES:
    d = w[w.name == s].copy()
    pre = int((d.post == 0).sum())
    if pre == 0:
        continue
    for y in ["diff", "logratio"]:
        b, se = contrast(d, y, "post")
        row = {"site": s, "outcome": y, "pre_months": pre, "post_months": int(d.post.sum()),
               "post_coef": b, "se_hac": se}
        # Placebo: pretend COD happened 12 months early, using pre-period data only.
        dp = d[d.post == 0].copy()
        if pre >= 24:
            dp["fake_post"] = (dp.rel_month >= -12).astype(int)
            pb, pse = contrast(dp, y, "fake_post")
            row.update(placebo_coef=pb, placebo_se=pse)
        rows.append(row)
res = pd.DataFrame(rows)
res.to_csv(NO2 / "site_prepost_contrasts.csv", index=False)
print("\n", res.round(3).to_string(index=False))


# 3. Figures --------------------------------------------------------------------------------
fig, axes = plt.subplots(4, 2, figsize=(10, 9.5), sharex=True)
for ax, s in zip(axes.flat, SITES):
    d = w[w.name == s]
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.plot(d.date, d["diff"], color=C2 if d.fuel.iloc[0] == "coal" else C1, lw=1.4)
    ax.plot(d.date, d["diff"].rolling(12, center=True).mean(), color=INK2, lw=1.2, ls="--")
    ax.axvline(d.commissioned.iloc[0], color=INK2, lw=1)
    a = inv.loc[s]
    ax.set_title(f"{s} ({a.fuel}, {a.mw} MW) — pre {a.pre_months} / post {a.post_months} mo",
                 fontsize=9)
for ax in axes[:, 0]:
    ax.set_ylabel("core − control\n(µmol/m²)")
fig.suptitle("Sentinel-5P NO₂, core minus control ring. Vertical line = commissioning; "
             "dashed = 12-month rolling mean", fontweight="bold", fontsize=10)
fig.tight_layout()
fig.savefig(NO2 / "fig1_core_minus_control_by_site.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 3.4))
seas = w.groupby(["fuel", "month"])[["core", "control"]].mean().reset_index()
for fuel, col in [("coal", C2), ("hydro", C1)]:
    sd = seas[seas.fuel == fuel]
    ax.plot(sd.month, sd.core, color=col, label=f"{fuel} core")
    ax.plot(sd.month, sd.control, color=col, ls=":", label=f"{fuel} control")
ax.set_xticks(range(1, 13))
ax.set_xlabel("Month")
ax.set_ylabel("NO₂ (µmol/m²)")
ax.set_title("Seasonal cycle: winter peak is larger than most plausible plant effects")
ax.legend(ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig(NO2 / "fig2_seasonality.png")
plt.close(fig)
print("\nFigures written to", NO2)
