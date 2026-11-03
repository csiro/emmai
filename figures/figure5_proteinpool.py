from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

WORKBOOK = Path("C:/Users/mol131/OneDrive - CSIRO/ecGEMs_Paper/Submission_files/Supplementary file 6.xlsx")
OUTPUT_PREFIX = Path("protein_pool_sensitivity")

POOL_FACTORS = np.array([0.50, 0.75, 1.00, 1.25, 1.50])

BASELINE_COLOR = "#6E6E6E"
EMMAI_COLOR = "#0072B2"
EMMAI_LIGHT = "#56B4E9"
GRID_COLOR = "#D9D9D9"

benchmarks = [
    {
        "label": "E. coli\niJO1366",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "columns": [f"iJO1366_pool_{factor:.2f}x" for factor in POOL_FACTORS],
        "quantitative": True,
        "marker": "o",
        "color": "#0072B2",
    },
    {
        "label": "E. coli\niML1515",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "columns": [f"iML1515_pool_{factor:.2f}x" for factor in POOL_FACTORS],
        "quantitative": True,
        "marker": "s",
        "color": "#56B4E9",
    },
    {
        "label": "P. aeruginosa\nreference",
        "sheet": "P. aeruginosa (reference)",
        "observed": "Biolog readings [2]",
        "columns": [f"pool_{factor:.2f}x" for factor in POOL_FACTORS],
        "quantitative": False,
        "marker": "o",
        "color": "#6E6E6E",
    },
    {
        "label": "P. aeruginosa\ndraft",
        "sheet": "P. aeruginosa (draft)",
        "observed": "Biolog readings [2]",
        "columns": [f"pool_{factor:.2f}x" for factor in POOL_FACTORS],
        "quantitative": False,
        "marker": "s",
        "color": "#9A9A9A",
    },
    {
        "label": "B. altitudinis\nBALT",
        "sheet": "B. altitudinis",
        "observed": "mean_max_rate",
        "columns": [f"pool_{factor:.2f}x" for factor in POOL_FACTORS],
        "quantitative": True,
        "marker": "o",
        "color": "#009E73",
    },
    {
        "label": "P. rettgeri\nPRET",
        "sheet": "P. rettgeri",
        "observed": "mean_max_rate",
        "columns": [f"pool_{factor:.2f}x" for factor in POOL_FACTORS],
        "quantitative": True,
        "marker": "s",
        "color": "#D55E00",
    },
]

sheets = pd.read_excel(WORKBOOK, sheet_name=None)
metric_rows = []

for benchmark in benchmarks:
    df = sheets[benchmark["sheet"]]
    observed = pd.to_numeric(df[benchmark["observed"]], errors="coerce").to_numpy(float)

    for factor, column in zip(POOL_FACTORS, benchmark["columns"]):
        predicted = pd.to_numeric(df[column], errors="coerce").to_numpy(float)
        mask = np.isfinite(observed) & np.isfinite(predicted)
        obs = observed[mask]
        pred = predicted[mask]

        pearson_r = pearsonr(obs, pred).statistic
        spearman_rho = spearmanr(obs, pred).statistic

        if benchmark["quantitative"]:
            residual = pred - obs
            mae = np.mean(np.abs(residual))
            rmse = np.sqrt(np.mean(residual ** 2))
        else:
            mae = np.nan
            rmse = np.nan

        metric_rows.append({
            "benchmark": benchmark["label"].replace("\n", " "),
            "pool_factor": factor, "pearson_r": pearson_r,
            "spearman_rho": spearman_rho, "mae": mae, "rmse": rmse
        })

metrics = pd.DataFrame(metric_rows)
metrics.to_csv(OUTPUT_PREFIX.with_name(OUTPUT_PREFIX.name + "_metrics.csv"), index=False)


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.7, alpha=0.65, zorder=0)
    ax.tick_params(labelsize=8)
    ax.axvline(1.0, color="#777777", lw=0.9, ls="--", zorder=1)


fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2))
plt.subplots_adjust(left=0.065, right=0.82, bottom=0.20, top=0.94, wspace=0.34)

ax = axes[0]

for benchmark in benchmarks:
    label = benchmark["label"].replace("\n", " ")
    sub = metrics[metrics["benchmark"] == label]

    ax.plot(
        sub["pool_factor"], sub["spearman_rho"], marker=benchmark["marker"],
        ms=5.2, lw=1.3, color=benchmark["color"], markeredgewidth=0,
        label=benchmark["label"], zorder=3,
    )

ax.set_xlabel("Protein pool multiplier", fontsize=9)
ax.set_ylabel(r"Spearman $\rho$", fontsize=9)
ax.set_xticks(POOL_FACTORS)
ax.set_ylim(-0.2, 0.75)
style_axis(ax)
shared_handles, shared_labels = ax.get_legend_handles_labels()
ax.text(-0.13, 1.08, "A", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")

ax = axes[1]

for benchmark in benchmarks:
    if not benchmark["quantitative"]:
        continue

    label = benchmark["label"].replace("\n", " ")
    sub = metrics[metrics["benchmark"] == label]

    ax.plot(sub["pool_factor"], sub["mae"], marker=benchmark["marker"], ms=5.2,
        lw=1.3, color=benchmark["color"], markeredgewidth=0, label=benchmark["label"],
        zorder=3)

ax.set_xlabel("Protein pool multiplier", fontsize=9)
ax.set_ylabel(r"MAE (h$^{-1}$)", fontsize=9)
ax.set_xticks(POOL_FACTORS)
style_axis(ax)
ax.text(-0.13, 1.08, "B", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")


ax = axes[2]

for benchmark in benchmarks:
    if not benchmark["quantitative"]:
        continue

    label = benchmark["label"].replace("\n", " ")
    sub = metrics[metrics["benchmark"] == label]

    ax.plot( sub["pool_factor"], sub["rmse"], marker=benchmark["marker"],
        ms=5.2, lw=1.3, color=benchmark["color"], markeredgewidth=0,
        label=benchmark["label"], zorder=3)

ax.set_xlabel("Protein pool multiplier", fontsize=9)
ax.set_ylabel(r"RMSE (h$^{-1}$)", fontsize=9)
ax.set_xticks(POOL_FACTORS)
style_axis(ax)
ax.text(-0.13, 1.08, "C", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")

fig.legend(
    shared_handles, shared_labels, loc="center left", bbox_to_anchor=(0.835, 0.56),
    frameon=False, fontsize=7.4, handlelength=2.2, handletextpad=0.6, labelspacing=0.75
)

fig.savefig(OUTPUT_PREFIX.with_suffix(".png"), dpi=400, bbox_inches="tight")
fig.savefig(OUTPUT_PREFIX.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(OUTPUT_PREFIX.with_suffix(".svg"), bbox_inches="tight")
plt.show()
