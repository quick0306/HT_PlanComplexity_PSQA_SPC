"""Generate summary figures for the HT PSQA analysis.

Outputs (saved as vector PDF and 300-DPI PNG):
- complexity_heatmap: Complexity metric correlation heatmap with cluster outlines
- site_distributions: Site-specific and pooled gamma index distributions (boxplots + jitter)
- ttdf_associations: Paired TTDF associations and global-local correlation differences
- site_associations: Within-site TTDF associations
- criterion_sensitivity: Gamma criterion sensitivity forest plots
"""

from argparse import ArgumentParser
from itertools import product
import os
from pathlib import Path
import tempfile

# Configure Matplotlib cache directory in the system temp folder to keep the workspace clean
ROOT = Path(__file__).resolve().parents[1]
cache_dir = Path(tempfile.gettempdir()) / "bpex_matplotlib_cache"
cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


SITES = ["Brain", "Head & neck", "Thorax", "Breast", "Abdomen", "Pelvis"]
CRITERIA = ["3%/3mm", "3%/2mm", "2%/2mm", "2%/1mm"]
PRIMARY = CRITERIA[1:3]
METRICS = {"GPR": ("GPR (%)", "pass%"), "mean_gamma": ("Mean γ", "mean_gamma")}
SERIES = [("Global", "#245789", "o"), ("Local", "#AA5000", "s")]


def get_gamma_col(norm: str, criterion: str, metric: str) -> str:
    """Construct column name for gamma metric."""
    return f"{norm.lower()}_{criterion.replace('/', '_')}_{METRICS[metric][1]}"


def format_criterion(criterion: str) -> str:
    """Format criterion string for display (e.g. '3%/2 mm')."""
    return criterion.replace("mm", " mm")


def save_figure(fig: plt.Figure, out_dir: Path, name: str) -> None:
    """Export figure in both vector PDF and 300-DPI raster PNG formats."""
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{name}.{ext}", dpi=300)
    plt.close(fig)


# ==============================================================================
# Complexity Heatmap & Complete-Linkage Clusters
# ==============================================================================
def plot_complexity_heatmap(results_dir: Path, out_dir: Path) -> None:
    """Heatmap of Spearman rho across 10 complexity metrics with cluster boxes."""
    corr = pd.read_excel(results_dir / "spearman_matrix.xlsx", index_col=0)
    clusters = pd.read_excel(results_dir / "clusters.xlsx").set_index("metric")["cluster"]

    fig, ax = plt.subplots(figsize=(7.2, 6.1), layout="constrained")
    edges = np.arange(len(corr) + 1) - 0.5
    palette = ax.pcolormesh(edges, edges, corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_aspect("equal")
    ax.invert_yaxis()

    ticks = np.arange(len(corr))
    ax.set_xticks(ticks, corr.columns, rotation=40, ha="right")
    ax.set_yticks(ticks, corr.index)

    # Annotate correlation values in cells
    for (row, col), val in np.ndenumerate(corr.to_numpy()):
        text_color = "white" if abs(val) >= 0.65 else "#202020"
        ax.text(col, row, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=8)

    # Draw outline boxes around complete-linkage clusters (cut height 0.50)
    cluster_order = clusters.reindex(corr.index).to_numpy()
    for cluster_id in np.unique(cluster_order):
        members = np.flatnonzero(cluster_order == cluster_id)
        if len(members) > 1:
            start = members[0] - 0.5
            ax.add_patch(Rectangle((start, start), len(members), len(members),
                                   fill=False, edgecolor="#202020", linewidth=1.5))

    fig.colorbar(palette, ax=ax, label="Spearman ρ", shrink=0.76, pad=0.025)
    save_figure(fig, out_dir, "complexity_heatmap")


# ==============================================================================
# Gamma Distributions by Treatment Site
# ==============================================================================
def plot_site_distributions(cohort: pd.DataFrame, results_dir: Path, out_dir: Path) -> None:
    """Boxplots and jittered points of GPR and Mean Gamma across 6 anatomical sites."""
    summary = pd.read_excel(results_dir / "gamma_distributions.xlsx")
    pooled = summary[summary["site"] == "All (pooled)"].set_index("metric")
    rng = np.random.default_rng(20261116)
    site_positions = np.arange(len(SITES))

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.7), layout="constrained")
    for panel, ((metric_key, criterion), ax) in enumerate(zip(product(METRICS, PRIMARY), axes.flat)):
        legend_labels = []
        for offset, (norm, color, marker) in zip((-0.17, 0.17), SERIES):
            col = get_gamma_col(norm, criterion, metric_key)
            values = [cohort.loc[cohort["site_en"] == s, col].to_numpy() for s in SITES]
            positions = site_positions + offset

            # Jittered scatter plot of individual plans
            for pos, group in zip(positions, values):
                jitter = rng.uniform(-0.052, 0.052, len(group))
                ax.scatter(pos + jitter, group, s=7, c=color, alpha=0.31, linewidths=0, zorder=2)

            # Boxplot showing median and IQR
            ax.boxplot(values, positions=positions, widths=0.25, showfliers=False,
                       boxprops={"color": color}, capprops={"color": color},
                       whiskerprops={"color": color}, medianprops={"color": color, "linewidth": 1.5},
                       zorder=3, manage_ticks=False)
            ax.scatter([], [], c=color, marker=marker, s=27, label=norm)

            # Extract pooled cohort skewness and ceiling rate
            item = pooled.loc[col]
            skew_val = item.get("skewness_uncorrected", item["skewness"])
            lbl = f"{norm}: skew {skew_val:.2f}"
            if metric_key == "GPR":
                lbl += f", ≥99.5%: {item['ceiling_pct']:.1f}%"
            legend_labels.append(lbl)

        if metric_key == "GPR":
            ax.axhline(99.5, color="#6C6C6C", linestyle="--", linewidth=0.75, zorder=1)

        ax.set_xticks(site_positions, SITES, rotation=20, ha="right")
        ax.set_ylabel(METRICS[metric_key][0])
        ax.set_title(f"{chr(65 + panel)}    {format_criterion(criterion)}\n" + "    |    ".join(legend_labels),
                     fontsize=9, loc="left")
        ax.grid(axis="y", color="#E1E1E1", linewidth=0.5, zorder=0)

    axes[0, 0].legend(loc="lower right", frameon=False)
    save_figure(fig, out_dir, "site_distributions")


# ==============================================================================
# Paired TTDF Associations (Global vs. Local)
# ==============================================================================
def draw_forest_intervals(ax: plt.Axes, positions: np.ndarray, estimates: np.ndarray,
                          low: np.ndarray, high: np.ndarray, color: str, marker: str,
                          label: str, *, vertical: bool = False) -> None:
    """Helper to draw forest plot point estimates and error bars."""
    est = np.asarray(estimates)
    err = [est - np.asarray(low), np.asarray(high) - est]
    kwargs = dict(color=color, fmt=marker, markersize=5, linewidth=1.15, capsize=2.4, label=label)
    if vertical:
        ax.errorbar(positions, est, yerr=err, **kwargs)
    else:
        ax.errorbar(est, positions, xerr=err, **kwargs)


def format_forest_axis(ax: plt.Axes, positions: np.ndarray, labels: list[str], xlabel: str) -> None:
    """Format vertical axis labels and zero reference line for forest plots."""
    ax.set_yticks(positions, labels)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.axvline(0, linestyle="--", color="#727272", linewidth=0.8, zorder=0)
    ax.grid(axis="y", color="#E4E4E4", linewidth=0.6, zorder=0)
    ax.legend(loc="best", frameon=False)


def plot_ttdf_associations(paired_df: pd.DataFrame, out_dir: Path) -> None:
    """Forest plots of TTDF associations by normalization and paired difference Δρ."""
    keys = [(metric, criterion) for metric in METRICS for criterion in PRIMARY]
    labels = [f"{METRICS[m][0].replace(' (%)', '')}, {format_criterion(c)}" for m, c in keys]
    positions = np.arange(len(keys))

    def get_subset(adj: str) -> pd.DataFrame:
        return paired_df[paired_df["adjustment"] == adj].set_index(["metric", "criterion"]).loc[keys]

    unadj = get_subset("unadjusted")
    site_adj = get_subset("site_adjusted")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.6, 4.2), layout="constrained")

    # Panel A: Normalization comparison (Global vs Local)
    for offset, (norm, color, marker) in zip((-0.11, 0.11), SERIES):
        col_base = f"rho_{norm.lower()}"
        draw_forest_intervals(ax_a, positions + offset, unadj[col_base],
                              unadj[f"{col_base}_ci_low"], unadj[f"{col_base}_ci_high"],
                              color, marker, norm)

    # Panel B: Paired differences delta_rho (global - local)
    for data, offset, color, marker, lbl in [
        (unadj, -0.11, "#454545", "o", "Unadjusted"),
        (site_adj, 0.11, "#177A82", "s", "Site-adjusted"),
    ]:
        draw_forest_intervals(ax_b, positions + offset, data["delta"],
                              data["delta_ci_low"], data["delta_ci_high"],
                              color, marker, lbl)

    format_forest_axis(ax_a, positions, labels, "Spearman ρ with TTDF")
    format_forest_axis(ax_b, positions, labels, "Δρ (global − local)")
    ax_a.set_title("A    Associations by normalization", loc="left")
    ax_b.set_title("B    Paired correlation differences", loc="left")
    save_figure(fig, out_dir, "ttdf_associations")


# ==============================================================================
# Within-Site TTDF Associations
# ==============================================================================
def plot_site_associations(results_dir: Path, out_dir: Path) -> None:
    """Forest plot of within-site Spearman rho with TTDF for 3%/2mm GPR."""
    site_df = pd.read_excel(results_dir / "site_associations.xlsx")
    positions = np.arange(len(SITES))
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")

    for offset, (norm, color, marker) in zip((-0.11, 0.11), SERIES):
        col_name = get_gamma_col(norm, PRIMARY[0], "GPR")
        data = site_df[site_df["outcome"] == col_name].set_index("site").loc[SITES]
        draw_forest_intervals(ax, positions + offset, data["rho"], data["lo"], data["hi"],
                              color, marker, norm)

    labels = [f"{s} (n = {n})" for s, n in zip(SITES, data["n"])]
    format_forest_axis(ax, positions, labels, "Spearman ρ with TTDF")
    save_figure(fig, out_dir, "site_associations")


# ==============================================================================
# Criterion Sensitivity Forest Plot
# ==============================================================================
def plot_criterion_sensitivity(paired_df: pd.DataFrame, out_dir: Path) -> None:
    """Forest plot of Spearman rho across all 4 dose/distance gamma criteria."""
    data_all = paired_df[paired_df["adjustment"] == "unadjusted"].set_index(["metric", "criterion"])
    positions = np.arange(len(CRITERIA))

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.1), layout="constrained", sharey=True)
    for panel, (metric_key, ax) in enumerate(zip(METRICS, axes)):
        subset = data_all.loc[[(metric_key, c) for c in CRITERIA]]
        for offset, (norm, color, marker) in zip((-0.045, 0.045), SERIES):
            base_col = f"rho_{norm.lower()}"
            draw_forest_intervals(ax, positions + offset, subset[base_col],
                                  subset[f"{base_col}_ci_low"], subset[f"{base_col}_ci_high"],
                                  color, marker, norm, vertical=True)

        ax.axhline(0, linestyle="--", color="#727272", linewidth=0.8)
        ax.set_xticks(positions, [format_criterion(c) for c in CRITERIA])
        ax.set_title(f"{chr(65 + panel)}    {METRICS[metric_key][0].replace(' (%)', '')}", loc="left")
        ax.set_xlabel("Gamma criterion")
        ax.legend(loc="best", frameon=False)

    axes[0].set_ylabel("Spearman ρ with TTDF")
    save_figure(fig, out_dir, "criterion_sensitivity")


# ==============================================================================
# CLI Entrypoint
# ==============================================================================
def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "analysis_cohort.csv")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "figures")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 10,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
    })

    cohort = pd.read_csv(args.data, float_precision="round_trip").query("analysis_cohort == 'Primary standard-site'")
    paired = pd.read_excel(args.results / "ttdf_paired_and_adjusted.xlsx")

    plot_complexity_heatmap(args.results, args.out)
    plot_site_distributions(cohort, args.results, args.out)
    plot_ttdf_associations(paired, args.out)
    plot_site_associations(args.results, args.out)
    plot_criterion_sensitivity(paired, args.out)


if __name__ == "__main__":
    main()
