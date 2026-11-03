from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

BASELINE_COLOR = "#6E6E6E"
EMMAI_COLOR = "#0072B2"
EMMAI_LIGHT = "#56B4E9"
SHUFFLE_FILL = "#E6E6E6"
SHUFFLE_EDGE = "#8A8A8A"
ZERO_COLOR = "#555555"
GRID_COLOR = "#D9D9D9"

ROBUSTNESS = [
    {
        "key": "ecoli_iml1515",
        "label": "E. coli\niML1515",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "gem10": "iML1515",
        "ec10": "eciML1515_10",
        "ec1000": "eciML1515",
        "quantitative": True,
    },
    {
        "key": "ecoli_ijo1366",
        "label": "E. coli\niJO1366",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "gem10": "iJO1366",
        "ec10": "eciJO1366_10",
        "ec1000": "eciJO1366",
        "quantitative": True,
    },
    {
        "key": "pao1_reference",
        "label": "P. aeruginosa\niMO1056",
        "sheet": "P. aeruginosa (reference)",
        "observed": "Biolog readings [2]",
        "gem10": "iMO1056",
        "ec10": "eciMO1056_10",
        "ec1000": "eciMO1056",
        "quantitative": False,
    },
    {
        "key": "pao1_draft",
        "label": "P. aeruginosa\nPAO1",
        "sheet": "P. aeruginosa (draft)",
        "observed": "Biolog readings [2]",
        "gem10": "PAO1",
        "ec10": "ecPAO1_10",
        "ec1000": "ecPAO1",
        "quantitative": False,
    },
    {
        "key": "balt",
        "label": "B. altitudinis\nBALT",
        "sheet": "B. altitudinis",
        "observed": "mean max rate",
        "gem10": "BALT",
        "ec10": "ecBALT_10",
        "ec1000": "ecBALT",
        "quantitative": True,
    },
    {
        "key": "pret",
        "label": "P. rettgeri\nPRET",
        "sheet": "P. rettgeri",
        "observed": "mean max rate",
        "gem10": "PRET",
        "ec10": "ecPRET_10",
        "ec1000": "ecPRET",
        "quantitative": True,
    },
]

ABLATION_BENCHMARK_MAP = {
    "E. coli iML1515": "E. coli\niML1515",
    "E. coli iJO1366": "E. coli\niJO1366",
    "P. aeruginosa iMO1056": "P. aeruginosa\niMO1056",
    "P. aeruginosa PAO1": "P. aeruginosa\nPAO1",
    "B. altitudinis": "B. altitudinis\nBALT",
    "P. rettgeri": "P. rettgeri\nPRET",
}

ABLATION_ORDER = list(ABLATION_BENCHMARK_MAP)
QUANT_ABLATION = ["E. coli iML1515", "E. coli iJO1366", "B. altitudinis", "P. rettgeri"]


ABLATION_SPECS = [
    {
        "benchmark": "E. coli iML1515",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "original": "emmai_iml",
        "uniform": "iml_uniform_median",
        "shuffle_prefix": "iml_shuffle_",
    },
    {
        "benchmark": "E. coli iJO1366",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "original": "emmai_ijo",
        "uniform": "ijo_uniform_median",
        "shuffle_prefix": "ijo_shuffle_",
    },
    {
        "benchmark": "P. aeruginosa iMO1056",
        "sheet": "P. aeruginosa (reference)",
        "observed": "Biolog readings [2]",
        "original": "emmai_imo",
        "uniform": "uniform_median",
        "shuffle_prefix": "shuffle_",
    },
    {
        "benchmark": "P. aeruginosa PAO1",
        "sheet": "P. aeruginosa (draft)",
        "observed": "Biolog readings [2]",
        "original": "emmai_pao1",
        "uniform": "uniform_median",
        "shuffle_prefix": "shuffle_",
    },
    {
        "benchmark": "B. altitudinis",
        "sheet": "B. altitudinis",
        "observed": "mean_max_rate",
        "original": "emmai_balt",
        "uniform": "uniform_median",
        "shuffle_prefix": "shuffle_",
    },
    {
        "benchmark": "P. rettgeri",
        "sheet": "P. rettgeri",
        "observed": "mean_max_rate",
        "original": "emmai_pret",
        "uniform": "uniform_median",
        "shuffle_prefix": "shuffle_",
    },
]

def clean_xy(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    return x[mask], y[mask]


def spearman_value(obs, pred):
    obs, pred = clean_xy(obs, pred)
    if len(obs) < 3:
        return np.nan
    return spearmanr(obs, pred).statistic


def bootstrap_spearman_ci(obs, pred, n_boot=10000, seed=20261002):
    obs, pred = clean_xy(obs, pred)
    n = len(obs)
    if n < 3:
        return np.nan, np.nan

    rng = np.random.default_rng(seed)
    vals = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            rho = spearmanr(obs[idx], pred[idx]).statistic
            if np.isfinite(rho):
                vals.append(rho)

    if not vals:
        return np.nan, np.nan
    return tuple(np.percentile(vals, [2.5, 97.5]))


def error_metrics(obs, pred):
    obs, pred = clean_xy(obs, pred)
    residual = pred - obs
    return {
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
    }


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.7, alpha=0.65, zorder=0)
    ax.tick_params(labelsize=8)


def add_panel_letter(ax, letter):
    ax.text(-0.13, 1.08, letter, transform=ax.transAxes,
            fontsize=13, fontweight="bold", fontstyle="italic", va="top", ha="left")


def resolve_sheet(sheets, spec):
    if "sheet" in spec:
        if spec["sheet"] not in sheets:
            raise KeyError(f"Sheet {spec['sheet']!r} not found. Available: {list(sheets)}")
        return spec["sheet"]
    

def load_robustness(workbook: Path):
    sheets = pd.read_excel(workbook, sheet_name=None)
    out = []

    for spec in ROBUSTNESS:
        sheet = resolve_sheet(sheets, spec)
        df = sheets[sheet].copy()
        required = [spec["observed"], spec["gem10"], spec["ec10"], spec["ec1000"]]
        
        sub = df[required].apply(pd.to_numeric, errors="coerce").dropna()
        item = dict(spec)
        item.update({
            "obs": sub[spec["observed"]].to_numpy(float),
            "gem10_values": sub[spec["gem10"]].to_numpy(float),
            "ec10_values": sub[spec["ec10"]].to_numpy(float),
            "ec1000_values": sub[spec["ec1000"]].to_numpy(float),
        })
        out.append(item)

    return out


def load_ablation(workbook):
    sheets = pd.read_excel(workbook, sheet_name=None)
    rows = []

    for spec in ABLATION_SPECS:
        df = sheets[spec["sheet"]].copy()
        shuffle_columns = sorted(
            c for c in df.columns if str(c).startswith(spec["shuffle_prefix"])
        )

        conditions = [
            ("original_EMMAi", spec["original"], np.nan),
            ("uniform_median", spec["uniform"], np.nan),
        ]

        for replicate, column in enumerate(shuffle_columns, start=1):
            conditions.append(("shuffle", column, replicate))

        for condition, column, replicate in conditions:
            obs, pred = clean_xy(
                pd.to_numeric(df[spec["observed"]], errors="coerce"),
                pd.to_numeric(df[column], errors="coerce"),
            )

            errors = error_metrics(obs, pred)

            rows.append({
                "benchmark": spec["benchmark"],
                "condition": condition,
                "replicate": replicate,
                "n": len(obs),
                "spearman_rho": spearman_value(obs, pred),
                "mae": errors["mae"],
                "rmse": errors["rmse"],
            })

    return pd.DataFrame(rows)


# figure
def plot_panel_a(ax, data, n_boot):
    x = np.arange(len(data), dtype=float)
    offsets = {"gem10": -0.16, "ec10": 0.0, "ec1000": 0.16}
    styles = {
        "gem10": (BASELINE_COLOR, "o", "GEM(10)"),
        "ec10": (EMMAI_LIGHT, "s", "ecGEM(10)"),
        "ec1000": (EMMAI_COLOR, "o", "ecGEM(1000)"),
    }

    for i, d in enumerate(data):
        xs = []
        ys = []

        for j, key in enumerate(["gem10", "ec10", "ec1000"]):
            pred = d[f"{key}_values"]
            rho = spearman_value(d["obs"], pred)
            lo, hi = bootstrap_spearman_ci(d["obs"], pred, n_boot=n_boot,
                                            seed=1000 + i * 10 + j)
            color, marker, _ = styles[key]
            xpos = x[i] + offsets[key]
            xs.append(xpos)
            ys.append(rho)
            ax.errorbar(
                xpos, rho,
                yerr=np.array([[rho - lo], [hi - rho]]),
                fmt=marker, ms=5.6, color=color,
                markeredgewidth=0, capsize=3.0, capthick=1.0, lw=1.0, zorder=3,
            )

        ax.plot(xs, ys, color="#B8B8B8", lw=1.2, zorder=1)

    ax.axhline(0, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([d["label"] for d in data], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Spearman $\rho$", fontsize=9)
    ax.set_ylim(-1.02, 1.02)
    style_axis(ax)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BASELINE_COLOR,
               markeredgewidth=0, label="GEM(10)", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=EMMAI_LIGHT,
               markeredgewidth=0, label="ecGEM(10)", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgewidth=0, label="ecGEM(1000)", markersize=5.5),
    ], frameon=False, fontsize=7.5, loc="lower right")


def plot_panel_b(ax, data):
    qdata = [d for d in data if d["quantitative"]]
    x = np.arange(len(qdata), dtype=float)
    cond_offsets = {"gem10": -0.18, "ec10": 0.0, "ec1000": 0.18}
    colors = {"gem10": BASELINE_COLOR, "ec10": EMMAI_LIGHT, "ec1000": EMMAI_COLOR}

    for i, d in enumerate(qdata):
        mae_x = []
        mae_y = []
        rmse_x = []
        rmse_y = []

        for key in ["gem10", "ec10", "ec1000"]:
            metrics = error_metrics(d["obs"], d[f"{key}_values"])

            mae_x.append(x[i] + cond_offsets[key] - 0.035)
            mae_y.append(metrics["mae"])
            rmse_x.append(x[i] + cond_offsets[key] + 0.035)
            rmse_y.append(metrics["rmse"])

            ax.scatter(mae_x[-1], mae_y[-1],
                       marker="o", s=29, color=colors[key], edgecolors="none", zorder=3)
            ax.scatter(rmse_x[-1], rmse_y[-1],
                       marker="s", s=29, color=colors[key], edgecolors="none", zorder=3)

        ax.plot(mae_x, mae_y, color="#B8B8B8", lw=1.2, zorder=1)
        ax.plot(rmse_x, rmse_y, color="#B8B8B8", lw=1.2, zorder=1)

    ax.set_xticks(x)
    ax.set_xticklabels([d["label"] for d in qdata], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Absolute prediction error (h$^{-1}$)", fontsize=9)
    style_axis(ax)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777",
               markeredgewidth=0, label="MAE", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#777777",
               markeredgewidth=0, label="RMSE", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BASELINE_COLOR,
               markeredgewidth=0, label="GEM(10)", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=EMMAI_LIGHT,
               markeredgewidth=0, label="ecGEM(10)", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgewidth=0, label="ecGEM(1000)", markersize=5.5),
    ], frameon=False, fontsize=7.1, ncol=2, loc="best")


def plot_panel_c(ax, data, seed=20261002):
    rng = np.random.default_rng(seed)
    values = []
    labels = []

    for d in data:
        delta = np.abs(d["ec1000_values"] - d["ec10_values"])
        values.append(delta)
        labels.append(d["label"])

    ax.boxplot(
        values,
        positions=np.arange(len(values)), widths=0.52,
        showfliers=False, patch_artist=True,
        medianprops=dict(color="#222222", linewidth=1.2),
        boxprops=dict(facecolor=SHUFFLE_FILL, edgecolor="#777777", linewidth=0.9),
        whiskerprops=dict(color="#777777", linewidth=0.9),
        capprops=dict(color="#777777", linewidth=0.9),
    )

    for i, delta in enumerate(values):
        jitter = rng.uniform(-0.17, 0.17, size=len(delta))
        ax.scatter(np.full(len(delta), i) + jitter, delta,
                   s=15, color=EMMAI_COLOR, alpha=0.72,
                   edgecolors="none", zorder=3)
        identical = int(np.sum(np.isclose(delta, 0.0, atol=1e-12, rtol=1e-9)))
        total = len(delta)
        ymax = np.max(delta) if len(delta) else 0
        ymin = np.min(delta) if len(delta) else 0
        ytxt = ymax + 0.06 * max(1e-9, ymax - ymin)
        ax.text(i, ytxt, f"{identical}/{total}", ha="center", va="bottom", fontsize=7.3)

    ax.axhline(0, color=ZERO_COLOR, lw=0.9, ls="--")
    ax.set_xticks(np.arange(len(values)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Absolute change in predicted biomass flux (h$^{-1}$)", fontsize=9)
    style_axis(ax)
    

def plot_panel_d(ax, ablation):
    x = np.arange(len(ABLATION_ORDER), dtype=float)
    shuffle_values = [
        ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "shuffle")]["spearman_rho"].dropna().to_numpy()
        for b in ABLATION_ORDER
    ]

    ax.boxplot(
        shuffle_values, positions=x, widths=0.50,
        showfliers=False, patch_artist=True,
        medianprops=dict(color="#222222", linewidth=1.1),
        boxprops=dict(facecolor=SHUFFLE_FILL, edgecolor=SHUFFLE_EDGE, linewidth=0.8),
        whiskerprops=dict(color=SHUFFLE_EDGE, linewidth=0.8),
        capprops=dict(color=SHUFFLE_EDGE, linewidth=0.8),
    )

    for i, b in enumerate(ABLATION_ORDER):
        original = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "original_EMMAi")]["spearman_rho"].iloc[0]
        uniform = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "uniform_median")]["spearman_rho"].iloc[0]
        ax.scatter(x[i] - 0.10, original, s=30, marker="o", color=EMMAI_COLOR, edgecolors="none", zorder=4)
        ax.scatter(x[i] + 0.10, uniform, s=30, marker="s", color=BASELINE_COLOR, edgecolors="none", zorder=4)

    ax.axhline(0, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([ABLATION_BENCHMARK_MAP[b] for b in ABLATION_ORDER], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Spearman $\rho$", fontsize=9)
    style_axis(ax)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgewidth=0, label="Original EMMAi", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=BASELINE_COLOR,
               markeredgewidth=0, label="Uniform median", markersize=5.5),
    ], frameon=False, fontsize=7.4, loc="best")


def plot_panel_e(ax, ablation):
    x = np.arange(len(QUANT_ABLATION), dtype=float)
    metric_shift = {"mae": -0.16, "rmse": 0.16}
    marker = {"mae": "o", "rmse": "s"}

    for metric in ["mae", "rmse"]:
        vals = [
            ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "shuffle")][metric].dropna().to_numpy()
            for b in QUANT_ABLATION
        ]
        ax.boxplot(
            vals, positions=x + metric_shift[metric], widths=0.22,
            showfliers=False, patch_artist=True,
            medianprops=dict(color="#222222", linewidth=1.0),
            boxprops=dict(facecolor=SHUFFLE_FILL, edgecolor=SHUFFLE_EDGE, linewidth=0.8),
            whiskerprops=dict(color=SHUFFLE_EDGE, linewidth=0.8),
            capprops=dict(color=SHUFFLE_EDGE, linewidth=0.8),
        )

        for i, b in enumerate(QUANT_ABLATION):
            original = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "original_EMMAi")][metric].iloc[0]
            uniform = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "uniform_median")][metric].iloc[0]
            ax.scatter(x[i] + metric_shift[metric] - 0.055, original,
                       s=29, marker=marker[metric], color=EMMAI_COLOR, edgecolors="none", zorder=4)
            ax.scatter(x[i] + metric_shift[metric] + 0.055, uniform,
                       s=29, marker=marker[metric], color=BASELINE_COLOR, edgecolors="none", zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels([ABLATION_BENCHMARK_MAP[b] for b in QUANT_ABLATION], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Absolute prediction error (h$^{-1}$)", fontsize=9)
    style_axis(ax)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777",
               markeredgewidth=0, label="MAE", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#777777",
               markeredgewidth=0, label="RMSE", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgewidth=0, label="Original EMMAi", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BASELINE_COLOR,
               markeredgewidth=0, label="Uniform median", markersize=5.5),
        Line2D([0], [0], color=SHUFFLE_EDGE, lw=6, alpha=0.35, label="50 shuffled models"),
    ], frameon=False, fontsize=7.0, ncol=2, loc="best")

def performance_percentile(original, shuffles, higher_is_better):
    shuffles = np.asarray(shuffles, dtype=float)
    shuffles = shuffles[np.isfinite(shuffles)]
    if len(shuffles) == 0 or not np.isfinite(original):
        return np.nan
    if higher_is_better:
        return 100.0 * np.mean(shuffles <= original)
    return 100.0 * np.mean(shuffles >= original)


def plot_panel_f(ax, ablation):
    x = np.arange(len(ABLATION_ORDER), dtype=float)
    spearman_p = []
    rmse_p = []

    for b in ABLATION_ORDER:
        original = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "original_EMMAi")].iloc[0]
        shuf = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "shuffle")]
        spearman_p.append(performance_percentile(original["spearman_rho"], shuf["spearman_rho"], True))
        if b in QUANT_ABLATION:
            rmse_p.append(performance_percentile(original["rmse"], shuf["rmse"], False))
        else:
            rmse_p.append(np.nan)

    ax.scatter(x - 0.07, spearman_p, marker="o", s=31, color=EMMAI_COLOR, edgecolors="none", zorder=3)
    valid = np.isfinite(rmse_p)
    ax.scatter(x[valid] + 0.07, np.asarray(rmse_p)[valid], marker="s", s=31,
               color=EMMAI_COLOR, edgecolors="none", zorder=3)

    ax.axhline(50, color="#777777", lw=0.9, ls="--")
    ax.set_ylim(-3, 103)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(x)
    ax.set_xticklabels([ABLATION_BENCHMARK_MAP[b] for b in ABLATION_ORDER], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("EMMAi performance percentile\nrelative to shuffled models", fontsize=9)
    style_axis(ax)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgewidth=0, label=r"Spearman $\rho$", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgewidth=0, label="RMSE", markersize=5.5),
    ], frameon=False, fontsize=7.5, loc="best")


def build_statistics(data, ablation, n_boot):
    rows = []

    for i, d in enumerate(data):
        for j, key in enumerate(["gem10", "ec10", "ec1000"]):
            pred = d[f"{key}_values"]
            rho = spearman_value(d["obs"], pred)
            lo, hi = bootstrap_spearman_ci(
                d["obs"], pred, n_boot=n_boot, seed=1000 + i * 10 + j
            )
            metrics = error_metrics(d["obs"], pred) if d["quantitative"] else {"mae": np.nan, "rmse": np.nan}
            rows.append({
                "analysis": "uptake_robustness",
                "benchmark": d["label"].replace("\n", " "),
                "condition": {"gem10": "GEM(10)", "ec10": "ecGEM(10)", "ec1000": "ecGEM(1000)"}[key],
                "replicate": np.nan,
                "n": len(d["obs"]),
                "spearman_rho": rho,
                "spearman_ci_low": lo,
                "spearman_ci_high": hi,
                "mae": metrics["mae"],
                "rmse": metrics["rmse"],
            })

        delta = np.abs(d["ec1000_values"] - d["ec10_values"])
        rows.append({
            "analysis": "uptake_sensitivity",
            "benchmark": d["label"].replace("\n", " "),
            "condition": "abs(ecGEM1000 - ecGEM10)",
            "replicate": np.nan,
            "n": len(delta),
            "mean_abs_change": float(np.mean(delta)),
            "median_abs_change": float(np.median(delta)),
            "max_abs_change": float(np.max(delta)),
            "unchanged_n": int(np.sum(np.isclose(delta, 0.0, atol=1e-12, rtol=1e-9))),
            "unchanged_total": len(delta),
        })

    for _, r in ablation.iterrows():
        rows.append({
            "analysis": "kcat_ablation",
            "benchmark": r["benchmark"],
            "condition": r["condition"],
            "replicate": r["replicate"],
            "n": r["n"],
            "spearman_rho": r["spearman_rho"],
            "mae": r["mae"],
            "rmse": r["rmse"],
        })

    for b in ABLATION_ORDER:
        original = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "original_EMMAi")].iloc[0]
        shuf = ablation[(ablation["benchmark"] == b) & (ablation["condition"] == "shuffle")]
        rows.append({
            "analysis": "kcat_ablation_percentile",
            "benchmark": b,
            "condition": "original_vs_shuffles",
            "replicate": np.nan,
            "n": original["n"],
            "spearman_percentile": performance_percentile(original["spearman_rho"], shuf["spearman_rho"], True),
            "rmse_percentile": (
                performance_percentile(original["rmse"], shuf["rmse"], False)
                if b in QUANT_ABLATION else np.nan
            ),
            "n_shuffles": len(shuf),
        })

    return pd.DataFrame(rows)

SUPP4_FILE = "Supplementary file 4.xlsx"
SUPP5_FILE = "Supplementary file 5.xlsx"
OUT_PREFIX = "emmai_figure4"
N_BOOT = 10000
DPI = 400

robustness = load_robustness(SUPP4_FILE)
ablation = load_ablation(SUPP5_FILE)

missing_ablation = [b for b in ABLATION_ORDER if b not in set(ablation["benchmark"])]
if missing_ablation:
    raise ValueError(
        f"Missing ablation benchmark(s): {missing_ablation}. "
        f"Loaded benchmark names: {sorted(ablation['benchmark'].unique())}"
    )

stats = build_statistics(robustness, ablation, N_BOOT)
stats.to_csv(f"{OUT_PREFIX}_statistics.csv", index=False)

fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.5))
plt.subplots_adjust(
    left=0.075, right=0.985, bottom=0.10, top=0.965,
    wspace=0.34, hspace=0.50,
)

plot_panel_a(axes[0, 0], robustness, N_BOOT)
plot_panel_b(axes[0, 1], robustness)
plot_panel_c(axes[0, 2], robustness)
plot_panel_d(axes[1, 0], ablation)
plot_panel_e(axes[1, 1], ablation)
plot_panel_f(axes[1, 2], ablation)

for letter, ax in zip("ABCDEF", axes.flat):
    add_panel_letter(ax, letter)

fig.savefig(f"{OUT_PREFIX}.png", dpi=DPI, bbox_inches="tight")
fig.savefig(f"{OUT_PREFIX}.pdf", bbox_inches="tight")
fig.savefig(f"{OUT_PREFIX}.svg", bbox_inches="tight")

plt.show()
