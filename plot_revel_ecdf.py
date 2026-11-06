"""
Compare REVEL score density distributions between:
  - Pathogenic variants filtered for any SpliceAI score > 0.5
  - VUS variants (already SpliceAI-filtered)

Used for Figure 4B (REVEL distribution, Pathogenic vs VUS).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import gaussian_kde, mannwhitneyu
import matplotlib as mpl
import os

# ── font ──────────────────────────────────────────────────────────────────────
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42

# ── paths ──────────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
PATH_PATHO = os.path.join(BASE, "input", "clinvar_pathogenic_missense_spliceai_revel.tsv")
PATH_VUS   = os.path.join(BASE, "input", "clinvar_spiceai_filtered05_revel.tsv")
OUT_PLOT   = os.path.join(BASE, "output", "revel_density_pathogenic_vs_vus.png")
OUT_STATS  = os.path.join(BASE, "output", "revel_stats_pathogenic_vs_vus.tsv")
os.makedirs(os.path.join(BASE, "output"), exist_ok=True)
REVEL_THR_PP3 = 0.644  # ClinGen SVI PP3 recommended cutoff
REVEL_THR_BP4 = 0.29   # ClinGen SVI BP4 recommended cutoff

# ── load & filter pathogenic ───────────────────────────────────────────────────
patho = pd.read_csv(PATH_PATHO, sep="\t", dtype=str)

# Numeric SpliceAI score columns already split
splice_cols = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]
for c in splice_cols:
    patho[c] = pd.to_numeric(patho[c], errors="coerce")

patho["REVEL"] = pd.to_numeric(patho["REVEL"], errors="coerce")

# Filter: any SpliceAI score > 0.5
mask = patho[splice_cols].gt(0.5).any(axis=1)
patho_filtered = patho[mask].copy()
print(f"Pathogenic variants total      : {len(patho)}")
print(f"Pathogenic after SpliceAI>0.5  : {len(patho_filtered)}")

revel_patho = patho_filtered["REVEL"].dropna()
print(f"  with REVEL score             : {len(revel_patho)}")

# ── load VUS ───────────────────────────────────────────────────────────────────
vus = pd.read_csv(PATH_VUS, sep="\t", dtype=str)
vus["REVEL"] = pd.to_numeric(vus["REVEL"], errors="coerce")
revel_vus = vus["REVEL"].dropna()
print(f"\nVUS variants total             : {len(vus)}")
print(f"  with REVEL score             : {len(revel_vus)}")

# ── export stats ─────────────────────────────────────────────────────────────
def group_stats(name, series, thr):
    return {
        "group":            name,
        "n":                len(series),
        "mean":             round(series.mean(), 4),
        "median":           round(series.median(), 4),
        "sd":               round(series.std(), 4),
        "pct25":            round(series.quantile(0.25), 4),
        "pct75":            round(series.quantile(0.75), 4),
        "n_above_thr":      int((series > thr).sum()),
        "pct_above_thr":    round(100 * (series > thr).mean(), 2),
    }

stat_rows = [
    group_stats("Pathogenic_SpliceAI>0.5", revel_patho, REVEL_THR_PP3),
    group_stats("VUS_SpliceAI>0.5",        revel_vus,   REVEL_THR_PP3),
]

# Mann-Whitney U test
U, pval = mannwhitneyu(revel_patho, revel_vus, alternative="two-sided")
print(f"\nMann-Whitney U = {U:.1f}, p = {pval:.3e}")

stats_df = pd.DataFrame(stat_rows)
stats_df["REVEL_threshold"] = REVEL_THR_PP3
stats_df["mannwhitney_U"]   = round(U, 1)
stats_df["mannwhitney_p"]   = f"{pval:.3e}"
stats_df.to_csv(OUT_STATS, sep="\t", index=False)
print(f"Stats saved to : {OUT_STATS}")
print(stats_df.to_string(index=False))

# ── compute KDE ───────────────────────────────────────────────────────────────
x_grid = np.linspace(0, 1, 500)

kde_p = gaussian_kde(revel_patho, bw_method=0.12)
kde_v = gaussian_kde(revel_vus,   bw_method=0.12)

# ── plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6))

ax.fill_between(x_grid, kde_p(x_grid), alpha=0.25, color="#d62728")
ax.plot(x_grid, kde_p(x_grid), color="#d62728", linewidth=2.2,
        label=f"Pathogenic (SpliceAI>0.5, n={len(revel_patho):,})")

ax.fill_between(x_grid, kde_v(x_grid), alpha=0.25, color="#1f77b4")
ax.plot(x_grid, kde_v(x_grid), color="#1f77b4", linewidth=2.2,
        label=f"VUS (SpliceAI>0.5, n={len(revel_vus):,})")

ax.axvline(REVEL_THR_PP3, color="#2ca02c", linestyle="--", linewidth=1.5,
           label=f"REVEL = {REVEL_THR_PP3} (PP3 cutoff)")
ax.axvline(REVEL_THR_BP4, color="#ff7f0e", linestyle="--", linewidth=1.5,
           label=f"REVEL = {REVEL_THR_BP4} (BP4 cutoff)")

ax.set_xlabel("REVEL score", fontsize=14)
ax.set_ylabel("Density", fontsize=14)
ax.set_title("REVEL score distribution\n(Pathogenic vs VUS, SpliceAI > 0.5)", fontsize=15)
ax.set_xlim(0, 1)
ax.tick_params(axis="both", labelsize=12)
ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.05))
ax.grid(True, which="major", linestyle="-",  alpha=0.3)
ax.grid(True, which="minor", linestyle=":",  alpha=0.15)
ax.legend(fontsize=12, loc="upper left")

plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=150)
plt.show()
print(f"\nPlot saved to: {OUT_PLOT}")
