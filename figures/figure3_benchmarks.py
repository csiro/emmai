from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr, pearsonr, kendalltau, linregress


BENCHMARKS = [
    {
        "key": "ecoli_iml1515",
        "label": "E. coli\niML1515",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "baseline": "iML1515",
        "emmai": "eciML1515",
        "baseline_label": "iML1515",
        "emmai_label": "EMMAi eciML1515",
        "quantitative": True,
        "y_label": r"Observed growth rate (h$^{-1}$)",
        "x_label": r"Predicted biomass flux (h$^{-1}$)",
    },
    {
        "key": "ecoli_ijo1366",
        "label": "E. coli\nAutoPACMEN",
        "sheet": "E. coli",
        "observed": "growth_rates",
        "baseline": "AutoPACMEN",
        "emmai": "eciJO1366",
        "baseline_label": "AutoPACMEN",
        "emmai_label": "EMMAi eciJO1366",
        "quantitative": True,
        "y_label": r"Observed growth rate (h$^{-1}$)",
        "x_label": r"Predicted biomass flux (h$^{-1}$)",
    },
    {
        "key": "pao1_reference",
        "label": "P. aeruginosa\niMO1056",
        "sheet": "P. aeruginosa (reference)",
        "observed": "Biolog readings [2]",
        "baseline": "iMO1056",
        "emmai": "eciMO1056",
        "baseline_label": "iMO1056",
        "emmai_label": "EMMAi eciMO1056",
        "quantitative": False,
        "y_label": "Biolog phenotype (a.u.)",
        "x_label": "Predicted biomass flux (h$^{-1}$)",
    },
    {
        "key": "pao1_draft",
        "label": "P. aeruginosa\nPAO1",
        "sheet": "P. aeruginosa (draft)",
        "observed": "Biolog readings [2]",
        "baseline": "PAO1",
        "emmai": "ecPAO1",
        "baseline_label": "PAO1",
        "emmai_label": "EMMAi ecPAO1",
        "quantitative": False,
        "y_label": "Biolog phenotype (a.u.)",
        "x_label": "Predicted biomass flux (h$^{-1}$)",
    },
    {
        "key": "balt",
        "label": "B. altitudinis\nBALT",
        "sheet": "B. altitudinis",
        "observed": "mean max rate",
        "baseline": "BALT",
        "emmai": "ecBALT",
        "baseline_label": "BALT",
        "emmai_label": "EMMAi ecBALT",
        "quantitative": True,
        "y_label": r"Apparent specific growth rate (h$^{-1}$)",
        "x_label": r"Predicted biomass flux (h$^{-1}$)",
    },
    {
        "key": "pret",
        "label": "P. rettgeri\nPRET",
        "sheet": "P. rettgeri",
        "observed": "mean max rate",
        "baseline": "PRET",
        "emmai": "ecPRET",
        "baseline_label": "PRET",
        "emmai_label": "EMMAi ecPRET",
        "quantitative": True,
        "y_label": r"Apparent specific growth rate (h$^{-1}$)",
        "x_label": r"Predicted biomass flux (h$^{-1}$)",
    },
]

# Colour palette
BASELINE_COLOR = "#6E6E6E"
EMMAI_COLOR = "#0072B2"
ZERO_COLOR = "#555555"
GRID_COLOR = "#D9D9D9"


def clean_xy(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    return x[mask], y[mask]


def metric_set(obs, pred):
    obs, pred = clean_xy(obs, pred)
    n = len(obs)
    if n < 2:
        return {
            "n": n,
            "pearson_r": np.nan,
            "spearman_rho": np.nan,
            "kendall_tau": np.nan,
            "rmse": np.nan,
            "mae": np.nan,
            "r_squared": np.nan,
        }

    pear = pearsonr(obs, pred)
    spear = spearmanr(obs, pred).statistic
    kend = kendalltau(obs, pred).statistic
    residual = pred - obs
    rmse = np.sqrt(np.mean(residual ** 2))
    mae = np.mean(np.abs(residual))

    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "n": n,
        "pearson_r": pear.statistic,
        "pearson_p": pear.pvalue,
        "spearman_rho": spear,
        "kendall_tau": kend,
        "rmse": rmse,
        "mae": mae,
        "r_squared": r2,
    }

def bootstrap_spearman_ci(obs, pred, n_boot=10000, seed=20260924):
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



def bootstrap_error_intervals(obs, baseline, emmai, n_boot=10000, seed=20261023):
    obs = np.asarray(obs, dtype=float)
    baseline = np.asarray(baseline, dtype=float)
    emmai = np.asarray(emmai, dtype=float)
    mask = np.isfinite(obs) & np.isfinite(baseline) & np.isfinite(emmai)
    obs, baseline, emmai = obs[mask], baseline[mask], emmai[mask]

    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(obs), size=(n_boot, len(obs)))
    errors_b = baseline - obs
    errors_e = emmai - obs

    # A single matched index matrix is used for both model predictions.
    mae_b = np.mean(np.abs(errors_b[indices]), axis=1)
    mae_e = np.mean(np.abs(errors_e[indices]), axis=1)
    rmse_b = np.sqrt(np.mean(errors_b[indices] ** 2, axis=1))
    rmse_e = np.sqrt(np.mean(errors_e[indices] ** 2, axis=1))

    result = {}
    for metric, baseline_boot, emmai_boot in (
        ("mae", mae_b, mae_e), ("rmse", rmse_b, rmse_e)
    ):
        delta_boot = baseline_boot - emmai_boot  # Positive is better for EMMAi.
        for name, values in (
            ("baseline", baseline_boot), ("emmai", emmai_boot), ("delta", delta_boot)
        ):
            low, high = np.percentile(values, [2.5, 97.5])
            result[f"{name}_{metric}_ci95_low"] = float(low)
            result[f"{name}_{metric}_ci95_high"] = float(high)
        result[f"delta_{metric}"] = float(
            (np.mean(np.abs(errors_b)) - np.mean(np.abs(errors_e)))
            if metric == "mae" else
            (np.sqrt(np.mean(errors_b ** 2)) - np.sqrt(np.mean(errors_e ** 2)))
        )
        result[f"delta_{metric}_ci_excludes_zero"] = bool(
            result[f"delta_{metric}_ci95_low"] > 0 or
            result[f"delta_{metric}_ci95_high"] < 0
        )
    return result


def paired_bootstrap_correlation_test(obs, baseline, emmai, method="spearman",
                                      n_boot=10000, seed=20261007):
    obs = np.asarray(obs, dtype=float)
    baseline = np.asarray(baseline, dtype=float)
    emmai = np.asarray(emmai, dtype=float)

    mask = np.isfinite(obs) & np.isfinite(baseline) & np.isfinite(emmai)
    obs = obs[mask]
    baseline = baseline[mask]
    emmai = emmai[mask]

    n = len(obs)
    if n < 3:
        return {
            "baseline": np.nan, "emmai": np.nan, "delta": np.nan,
            "ci95_low": np.nan, "ci95_high": np.nan, "bootstrap_probability_delta_gt_0": np.nan,
            "ci_excludes_zero": False, "n_valid_bootstrap": 0,
        }

    if method == "spearman":
        corr = lambda x, y: spearmanr(x, y).statistic
    elif method == "pearson":
        corr = lambda x, y: pearsonr(x, y).statistic

    baseline_corr = corr(obs, baseline)
    emmai_corr = corr(obs, emmai)
    observed_delta = emmai_corr - baseline_corr

    rng = np.random.default_rng(seed)
    deltas = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for _ in range(n_boot):
            idx = rng.integers(0, n, n)

            b = corr(obs[idx], baseline[idx])
            e = corr(obs[idx], emmai[idx])

            if np.isfinite(b) and np.isfinite(e):
                deltas.append(e - b)

    if not deltas:
        return {
            "baseline": baseline_corr,
            "emmai": emmai_corr,
            "delta": observed_delta,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
            "bootstrap_probability_delta_gt_0": np.nan,
            "ci_excludes_zero": False,
            "n_valid_bootstrap": 0,
        }

    deltas = np.asarray(deltas)
    ci_low, ci_high = np.percentile(deltas, [2.5, 97.5])

    return {
        "baseline": baseline_corr,
        "emmai": emmai_corr,
        "delta": observed_delta,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "bootstrap_probability_delta_gt_0": float(np.mean(deltas > 0)),
        "ci_excludes_zero": bool(ci_low > 0 or ci_high < 0),
        "n_valid_bootstrap": len(deltas),
    }


def fit_line(x, y):
    x, y = clean_xy(x, y)
    if len(x) < 2 or np.allclose(x, x[0]):
        return None
    fit = linregress(x, y)
    xx = np.linspace(np.min(x), np.max(x), 100)
    yy = fit.intercept + fit.slope * xx
    return xx, yy


def load_benchmark_data(workbook: Path):
    sheets = pd.read_excel(workbook, sheet_name=None)
    out = []

    for b in BENCHMARKS:

        df = sheets[b["sheet"]].copy()
        required = [b["observed"], b["baseline"], b["emmai"]]

        sub = df[required].apply(pd.to_numeric, errors="coerce").dropna()
        obs = sub[b["observed"]].to_numpy(float)
        baseline = sub[b["baseline"]].to_numpy(float)
        emmai = sub[b["emmai"]].to_numpy(float)

        item = dict(b)
        item.update({"obs": obs, "baseline_values": baseline, "emmai_values": emmai})
        out.append(item)

    return out


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.7, alpha=0.65, zorder=0)
    ax.tick_params(labelsize=8)


def add_panel_letter(ax, letter):
    ax.text(-0.13, 1.08, letter, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="top", ha="left")


def plot_panel_a(ax, data, n_boot):
    x = np.arange(len(data), dtype=float)
    offset = 0.11

    for i, d in enumerate(data):
        b = metric_set(d["obs"], d["baseline_values"])["spearman_rho"]
        e = metric_set(d["obs"], d["emmai_values"])["spearman_rho"]

        b_ci = bootstrap_spearman_ci(
            d["obs"], d["baseline_values"], n_boot=n_boot, seed=500 + i
        )
        e_ci = bootstrap_spearman_ci(
            d["obs"], d["emmai_values"], n_boot=n_boot, seed=1500 + i
        )

        paired = paired_bootstrap_correlation_test(
            d["obs"], d["baseline_values"], d["emmai_values"], method="spearman",
            n_boot=n_boot, seed=2500 + i
        )

        ax.plot([x[i] - offset, x[i] + offset], [b, e],
            color="#B8B8B8", lw=1.4, zorder=1)

        ax.errorbar(x[i] - offset, b, yerr=np.array([[b - b_ci[0]], [b_ci[1] - b]]),
            fmt="o",  ms=5.5, color=BASELINE_COLOR, capsize=2.5, lw=1.0,
            zorder=3)

        ax.errorbar(
            x[i] + offset, e, yerr=np.array([[e - e_ci[0]], [e_ci[1] - e]]),
            fmt="o", ms=5.5, color=EMMAI_COLOR, capsize=2.5, lw=1.0,
            zorder=3)

        # Asterisk = paired bootstrap 95% CI for delta rho excludes zero.
        if paired["ci_excludes_zero"]:
            upper = max(b_ci[1], e_ci[1])
            y_star = min(0.96, upper + 0.07)
            ax.text(x[i], y_star, "*", ha="center", va="center", fontsize=11, fontweight="bold",
                color="#333333", zorder=4)

    ax.axhline(0, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(
        [d["label"] for d in data], rotation=35, ha="right", fontsize=7.5
    )
    ax.set_ylabel(r"Spearman $\rho$", fontsize=9)
    ax.set_ylim(-1.02, 1.02)
    style_axis(ax)

    ax.legend(
        handles=[
            Line2D(
                [0], [0], marker="o", color="none",
                markerfacecolor=BASELINE_COLOR, markeredgecolor=BASELINE_COLOR,
                label="Baseline", markersize=5.5,
            ),
            Line2D(
                [0], [0], marker="o", color="none",
                markerfacecolor=EMMAI_COLOR, markeredgecolor=EMMAI_COLOR,
                label="EMMAi", markersize=5.5,
            ),
        ],
        frameon=False, fontsize=7.2, loc="lower right",
    )


def plot_panel_b(ax, data, n_boot):
    qdata = [d for d in data if d["quantitative"]]
    x = np.arange(len(qdata), dtype=float)
    offset = 0.10
    ymax = 0.0

    for i, d in enumerate(qdata):
        bm = metric_set(d["obs"], d["baseline_values"])
        em = metric_set(d["obs"], d["emmai_values"])
        ci = bootstrap_error_intervals(
            d["obs"], d["baseline_values"], d["emmai_values"],
            n_boot=n_boot, seed=9300 + i
        )

        for metric, marker, shift in [("mae", "o", -0.08), ("rmse", "s", +0.08)]:
            x0 = x[i] + shift - offset
            x1 = x[i] + shift + offset
            y0, y1 = bm[metric], em[metric]
            lo0, hi0 = ci[f"baseline_{metric}_ci95_low"], ci[f"baseline_{metric}_ci95_high"]
            lo1, hi1 = ci[f"emmai_{metric}_ci95_low"], ci[f"emmai_{metric}_ci95_high"]

            ax.plot([x0, x1], [y0, y1], color="#B8B8B8", lw=1.3, zorder=1)
            for xpos, value, lo, hi, color in (
                (x0, y0, lo0, hi0, BASELINE_COLOR),
                (x1, y1, lo1, hi1, EMMAI_COLOR),
            ):
                ax.errorbar(
                    xpos, value,
                    yerr=[[max(0, value - lo)], [max(0, hi - value)]],
                    fmt=marker, markersize=5.2, color=color,
                    capsize=2.5, elinewidth=1.0, capthick=1.0, zorder=3,
                )
            pair_top = max(hi0, hi1)
            ymax = max(ymax, pair_top)
            if ci[f"delta_{metric}_ci_excludes_zero"]:
                star_y = pair_top + 0.045  
                ax.text((x0 + x1) / 2, star_y, "*", ha="center", va="bottom",
                        fontsize=11, fontweight="bold", color="#333333", zorder=4)
                ymax = max(ymax, star_y + 0.035)

    ax.set_xticks(x)
    ax.set_xticklabels([d["label"] for d in qdata], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Absolute prediction error (h$^{-1}$)", fontsize=9)
    ax.set_ylim(bottom=0, top=max(ymax * 1.08, 0.1))
    style_axis(ax)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777",
               markeredgecolor="#777777", label="MAE", markersize=5.5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#777777",
               markeredgecolor="#777777", label="RMSE", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BASELINE_COLOR,
               markeredgecolor=BASELINE_COLOR, label="Baseline", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=EMMAI_COLOR,
               markeredgecolor=EMMAI_COLOR, label="EMMAi", markersize=5.5),
    ], frameon=False, fontsize=6.5, ncol=3, loc="upper center",
       bbox_to_anchor=(0.32, 0.95))


def plot_panel_c(ax, data, seed=20260924):
    qdata = [d for d in data if d["quantitative"]]
    rng = np.random.default_rng(seed)
    values = []

    for d in qdata:
        delta = np.abs(d["baseline_values"] - d["obs"]) - np.abs(d["emmai_values"] - d["obs"])
        values.append(delta)

    bp = ax.boxplot(values, positions=np.arange(len(values)), widths=0.52,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="#222222", linewidth=1.2),
                    boxprops=dict(facecolor="#E6E6E6", edgecolor="#777777", linewidth=0.9),
                    whiskerprops=dict(color="#777777", linewidth=0.9),
                    capprops=dict(color="#777777", linewidth=0.9))

    for i, delta in enumerate(values):
        jitter = rng.uniform(-0.17, 0.17, size=len(delta))
        point_colors = np.where(delta > 0, EMMAI_COLOR, BASELINE_COLOR)
        ax.scatter(np.full(len(delta), i) + jitter, delta, s=15, c=point_colors,
                   alpha=0.72, edgecolors="none", zorder=3)
        improved = int(np.sum(delta > 0))
        total = len(delta)
        ymax = np.max(delta) if len(delta) else 0
        ymin = np.min(delta) if len(delta) else 0
        ytxt = ymax + 0.06 * max(1e-9, ymax-ymin)
        ax.text(i, ytxt, f"{improved}/{total}", ha="center", va="bottom", fontsize=7.5)

    ax.axhline(0, color=ZERO_COLOR, lw=0.9, ls="--")
    ax.set_xticks(np.arange(len(qdata)))
    ax.set_xticklabels([d["label"] for d in qdata], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"Reduction in absolute error (h$^{-1}$)", fontsize=9)
    style_axis(ax)
    ax.text(0.02, 0.03, "", transform=ax.transAxes, va="bottom", ha="left", fontsize=7.2,
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="none", alpha=0.72))


def plot_scatter_panel(ax, d):
    obs = d["obs"]
    base = d["baseline_values"]
    emmai = d["emmai_values"]

    bm = metric_set(obs, base)
    em = metric_set(obs, emmai)

    base_legend = rf"{d['baseline_label']} ($r$={bm['pearson_r']:.2f})"
    emmai_legend = rf"{d['emmai_label']} ($r$={em['pearson_r']:.2f})"

    ax.scatter(base, obs, s=22, alpha=0.78, color=BASELINE_COLOR,
               edgecolors="white", linewidths=0.35, label=base_legend, zorder=3)
    ax.scatter(emmai, obs, s=22, alpha=0.82, color=EMMAI_COLOR,
               edgecolors="white", linewidths=0.35, label=emmai_legend, zorder=3)

    bfit = fit_line(base, obs)
    efit = fit_line(emmai, obs)
    if bfit:
        ax.plot(*bfit, color=BASELINE_COLOR, lw=1.0, alpha=0.85)
    if efit:
        ax.plot(*efit, color=EMMAI_COLOR, lw=1.0, alpha=0.9)

    if d["quantitative"]:
        x_max = max(np.max(base), np.max(emmai))
        y_max = np.max(obs)

        x_min = min(np.min(base), np.min(emmai))
        y_min = np.min(obs)

        x_pad = 0.05 * (x_max - x_min if x_max > x_min else 1.0)
        y_pad = 0.05 * (y_max - y_min if y_max > y_min else 1.0)

        x_lo = 0.0
        x_hi = x_max + x_pad
        y_lo = 0.0
        y_hi = y_max + y_pad

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)

        line_hi = min(x_hi, y_hi)

        ax.plot([0.0, line_hi], [0.0, line_hi], ls="--",
            lw=0.9, color="#222222", alpha=0.65, label="1:1 agreement")

    ax.set_xlabel(d["x_label"], fontsize=8.2)
    ax.set_ylabel(d["y_label"], fontsize=8.2)
    style_axis(ax)

    ax.ticklabel_format(style="plain", axis="both", useOffset=False)
    ax.legend(loc="best",
              frameon=False, fontsize=6.8, handletextpad=0.4,
              borderaxespad=0.3, labelspacing=0.25)


# Figure settings and statistics

workbook = Path("Supplementary file 4.xlsx")
out_prefix = "emmai_figure3"
n_bootstrap = 10000
dpi = 400

data = load_benchmark_data(workbook)
statistics_rows = []

for i, d in enumerate(data):
    bm = metric_set(d["obs"], d["baseline_values"])
    em = metric_set(d["obs"], d["emmai_values"])

    b_spear_ci = bootstrap_spearman_ci(
        d["obs"], d["baseline_values"],
        n_boot=n_bootstrap, seed=9000 + i * 10
    )

    e_spear_ci = bootstrap_spearman_ci(
        d["obs"], d["emmai_values"],
        n_boot=n_bootstrap, seed=9001 + i * 10
    )

    paired_spearman = paired_bootstrap_correlation_test(
        d["obs"], d["baseline_values"], d["emmai_values"],
        method="spearman", n_boot=n_bootstrap, seed=9100 + i
    )

    paired_pearson = paired_bootstrap_correlation_test(
        d["obs"], d["baseline_values"], d["emmai_values"],
        method="pearson", n_boot=n_bootstrap, seed=9200 + i
    )

    row = {
        "benchmark": d["key"],
        "benchmark_label": d["label"],
        "baseline_model": d["baseline_label"],
        "emmai_model": d["emmai_label"],
        "quantitative_growth_rate": d["quantitative"],
        "n": bm["n"],

        "baseline_pearson_r": bm["pearson_r"],
        "baseline_pearson_p": bm["pearson_p"],
        "emmai_pearson_r": em["pearson_r"],
        "emmai_pearson_p": em["pearson_p"],
        "delta_pearson_r": paired_pearson["delta"],
        "delta_pearson_ci95_low": paired_pearson["ci95_low"],
        "delta_pearson_ci95_high": paired_pearson["ci95_high"],
        "delta_pearson_bootstrap_probability_gt_0":
            paired_pearson["bootstrap_probability_delta_gt_0"],
        "delta_pearson_ci_excludes_zero":
            paired_pearson["ci_excludes_zero"],

        "baseline_spearman_rho": bm["spearman_rho"],
        "baseline_spearman_ci95_low": b_spear_ci[0],
        "baseline_spearman_ci95_high": b_spear_ci[1],
        "emmai_spearman_rho": em["spearman_rho"],
        "emmai_spearman_ci95_low": e_spear_ci[0],
        "emmai_spearman_ci95_high": e_spear_ci[1],
        "delta_spearman_rho": paired_spearman["delta"],
        "delta_spearman_ci95_low": paired_spearman["ci95_low"],
        "delta_spearman_ci95_high": paired_spearman["ci95_high"],
        "delta_spearman_bootstrap_probability_gt_0":
            paired_spearman["bootstrap_probability_delta_gt_0"],
        "delta_spearman_ci_excludes_zero":
            paired_spearman["ci_excludes_zero"],

        "baseline_kendall_tau": bm["kendall_tau"],
        "emmai_kendall_tau": em["kendall_tau"],
        "baseline_r_squared_agreement": bm["r_squared"],
        "emmai_r_squared_agreement": em["r_squared"],

        "baseline_mae": bm["mae"] if d["quantitative"] else np.nan,
        "emmai_mae": em["mae"] if d["quantitative"] else np.nan,
        "baseline_rmse": bm["rmse"] if d["quantitative"] else np.nan,
        "emmai_rmse": em["rmse"] if d["quantitative"] else np.nan,
    }

    if d["quantitative"]:
        row.update(bootstrap_error_intervals(
                d["obs"], d["baseline_values"], d["emmai_values"], n_boot=n_bootstrap, seed=9300 + i))
        delta_error = (np.abs(d["baseline_values"] - d["obs"]) - np.abs(d["emmai_values"] - d["obs"]))

        row.update({
            "substrates_improved": int(np.sum(delta_error > 0)),
            "substrates_total": len(delta_error),
            "substrates_improved_percent":
                100 * np.mean(delta_error > 0),
            "mean_reduction_in_absolute_error":
                float(np.mean(delta_error)),
            "median_reduction_in_absolute_error":
                float(np.median(delta_error)),
        })

    else:
        row.update({
            "substrates_improved": np.nan,
            "substrates_total": np.nan,
            "substrates_improved_percent": np.nan,
            "mean_reduction_in_absolute_error": np.nan,
            "median_reduction_in_absolute_error": np.nan,
        })

    statistics_rows.append(row)

statistics_df = pd.DataFrame(statistics_rows)
statistics_path = Path(f"{out_prefix}_statistics.csv")
statistics_df.to_csv(statistics_path, index=False)

fig, axes = plt.subplots(3, 3, figsize=(13.2, 12.2))

plt.subplots_adjust(left=0.075, right=0.985, bottom=0.07, top=0.955,
    wspace=0.34, hspace=0.48)

plot_panel_a(axes[0, 0], data, n_bootstrap)
plot_panel_b(axes[0, 1], data, n_bootstrap)
plot_panel_c(axes[0, 2], data)

for ax, d in zip(axes.flat[3:], data):
    plot_scatter_panel(ax, d)

for letter, ax in zip("ABCDEFGHI", axes.flat):
    add_panel_letter(ax, letter)

output = Path(out_prefix)

png = output.with_suffix(".png")
pdf = output.with_suffix(".pdf")
svg = output.with_suffix(".svg")

fig.savefig(png, dpi=dpi, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(svg, bbox_inches="tight")
plt.show()

print(f"Saved: {png}")
print(f"Saved: {pdf}")
print(f"Saved: {svg}")
print(f"Saved: {statistics_path}")


