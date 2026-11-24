"""Primary statistical analysis and process control limits for HT PSQA.

Outputs:
- Complexity metric clustering (clusters.xlsx, spearman_matrix.xlsx, linkage.xlsx)
- Gamma index and GPR distributions (gamma_distributions.xlsx)
- Primary TTDF global vs. local paired associations (ttdf_paired_and_adjusted.xlsx)
- Secondary complexity metric associations (secondary_associations.xlsx)
- Site-stratified correlations & leave-one-out diagnostics (site_associations.xlsx)
- CSI sensitivity analysis (csi_sensitivity.xlsx)
- Planning-parameter partial-rank associations (planning_parameter_associations.xlsx)
- Statistical process control limits: Shewhart vs. WSD (spc_limits*.xlsx, spc_validation_*.xlsx)
- Chronological validation split (temporal_associations.xlsx, chronological_split.xlsx)
"""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform
from statsmodels.stats.diagnostic import normal_ad

if __package__:
    from .result_tables import save_table
    from .study_settings import BOOTSTRAP_DRAWS, chronological_split
else:
    from result_tables import save_table
    from study_settings import BOOTSTRAP_DRAWS, chronological_split


# Analysis cohort configuration
SITES = ["Brain", "Head & neck", "Thorax", "Breast", "Abdomen", "Pelvis"]
METRICS = ["CFNS75", "MF", "mFLOT", "L0NS", "CLS", "fDISC", "PSTV", "MI", "TTDF", "centroid"]
SECONDARY_METRICS = ["CFNS75", "MI", "L0NS", "centroid"]
PRIMARY_CRITERIA = ["3%_2mm", "2%_2mm"]
ALL_CRITERIA = [*PRIMARY_CRITERIA, "3%_3mm", "2%_1mm"]

# Outcome definitions
PAIRS = [(criterion, kind) for kind in ["pass%", "mean_gamma"] for criterion in ALL_CRITERIA]
OUTCOMES = [f"{norm}_{criterion}_{kind}" for criterion, kind in PAIRS for norm in ["global", "local"]]
PRIMARY_OUTCOMES = [col for col in OUTCOMES if any(c in col for c in PRIMARY_CRITERIA)]
GLOBAL_GPR = "global_3%_2mm_pass%"
LOCAL_GPR = "local_3%_2mm_pass%"

# Random seeds & bootstrap iterations
SEED_MAIN = 20260929
SEED_SITE = 20260930
SEED_TEMPORAL = 20261116
PLANNING_PARAMETERS = ["FieldWidth", "Pitch", "MF", "TargetLength"]
ADJUSTMENTS = ["unadjusted", "site_adjusted", "site_fw_pitch_mf", "site_fw_pitch_mf_length"]


def rank_and_center(values: np.ndarray) -> np.ndarray:
    """Compute column-wise fractional ranks centered around the mean."""
    ranked = stats.rankdata(values, axis=0, method="average")
    return ranked - ranked.mean(axis=0)


def fast_correlation(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Compute Pearson correlation matrix of pre-centered rank matrices (Spearman rho)."""
    scales = np.sqrt(np.sum(x**2, axis=0)[:, None] * np.sum(y**2, axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (x.T @ y) / scales


def get_patient_clusters(patient_ids: pd.Series) -> list[np.ndarray]:
    """Extract list of row index arrays belonging to each unique patient ID."""
    if patient_ids.isna().any():
        raise ValueError("Patient IDs are required for patient-cluster bootstrap.")
    codes, labels = pd.factorize(patient_ids, sort=False)
    return [np.flatnonzero(codes == i) for i in range(len(labels))]


def resample_clusters(clusters: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    """Resample patient clusters with replacement, preserving within-patient plan correlation."""
    indices = rng.integers(0, len(clusters), size=len(clusters))
    return np.concatenate([clusters[i] for i in indices])


def quantile_ci(draws: np.ndarray) -> np.ndarray:
    """Calculate 95% empirical bootstrap confidence intervals (2.5% and 97.5% quantiles)."""
    return np.nanquantile(draws, [0.025, 0.975], axis=0)


# ==============================================================================
# 1. Complexity Metric Hierarchical Clustering
# ==============================================================================
def run_clustering(primary: pd.DataFrame, output: Path) -> None:
    """Hierarchical complete-linkage clustering of 10 plan complexity metrics.

    Uses distance d = 1 - |rho|, cut at height 0.50 to guarantee pairwise |rho| >= 0.50
    within each multi-member cluster.
    """
    corr = primary[METRICS].corr(method="spearman")
    distances = 1 - corr.abs().to_numpy()
    np.fill_diagonal(distances, 0)

    # Complete linkage clustering
    tree = linkage(squareform(distances), method="complete")
    order = leaves_list(tree)

    # Save ordered correlation matrix, dendrogram linkage, and cluster membership
    save_table(corr.iloc[order, order], output, "spearman_matrix", index_label="metric")
    save_table(pd.DataFrame(tree, columns=["left", "right", "distance", "size"]), output, "linkage")
    save_table({
        "metric": METRICS,
        "cluster": fcluster(tree, t=0.5, criterion="distance"),
        "retained": [m in ["TTDF", *SECONDARY_METRICS] for m in METRICS],
    }, output, "clusters")


# ==============================================================================
# 2. Gamma Index and Passing Rate Distributions
# ==============================================================================
def run_distributions(primary: pd.DataFrame, output: Path) -> None:
    """Calculate descriptive statistics, ceiling rates, and normality tests for gamma outcomes."""
    rows = []
    for outcome in PRIMARY_OUTCOMES:
        site_groups = [(site, primary.loc[primary.site_en == site, outcome]) for site in SITES]
        kruskal = stats.kruskal(*(vals for _, vals in site_groups))

        # Overall pooled cohort and per-site distribution
        for site, vals in [("All (pooled)", primary[outcome]), *site_groups]:
            ceiling_count = int((vals >= 99.5).sum()) if outcome.endswith("pass%") else None
            ad_stat, ad_p = normal_ad(vals.to_numpy(float))
            rows.append({
                "metric": outcome,
                "site": site,
                "n": len(vals),
                "mean": vals.mean(),
                "sd": vals.std(ddof=1),
                "skewness": stats.skew(vals, bias=False),
                "skewness_uncorrected": stats.skew(vals, bias=True),
                "anderson_darling_A2": ad_stat,
                "anderson_darling_p": ad_p,
                "ceiling_n": ceiling_count,
                "ceiling_pct": None if ceiling_count is None else 100.0 * ceiling_count / len(vals),
                "kruskal_H": kruskal.statistic if site == "All (pooled)" else None,
                "kruskal_p": kruskal.pvalue if site == "All (pooled)" else None,
            })
    save_table(rows, output, "gamma_distributions")


# ==============================================================================
# 3. Primary & Secondary Associations with Bootstrap
# ==============================================================================
def estimate_ranks_and_adjusted(values: np.ndarray, site_indices: np.ndarray):
    """Compute unadjusted ranks and site-mean-centered ranks for partial correlation."""
    ranked = rank_and_center(values)
    residual = ranked.copy()
    for site in np.unique(site_indices):
        mask = site_indices == site
        residual[mask] -= ranked[mask].mean(axis=0)

    n_metrics = 1 + len(SECONDARY_METRICS)
    base_corr = fast_correlation(ranked[:, :n_metrics], ranked[:, n_metrics:])
    adj_corr = fast_correlation(residual[:, :1], residual[:, n_metrics:])[0]
    return base_corr, adj_corr


def estimate_planning_adjusted(values: np.ndarray, site_design: np.ndarray,
                               parameters: np.ndarray, n_adjusters: int) -> np.ndarray:
    """Residualize TTDF and outcomes on site indicators and ranked planning parameters."""
    ranked_values = rank_and_center(values)
    ranked_parameters = rank_and_center(parameters[:, :n_adjusters])
    design = np.column_stack([site_design, ranked_parameters])
    fitted = design @ np.linalg.lstsq(design, ranked_values, rcond=None)[0]
    residual = ranked_values - fitted
    return fast_correlation(residual[:, :1], residual[:, 1:])[0]


def build_paired_rows(primary: pd.DataFrame, point: np.ndarray, draws: np.ndarray,
                      adjustments: list[str], pairs: list[tuple[str, str]]) -> list[dict]:
    """Format paired global vs. local TTDF correlation comparison table."""
    rows = []
    for pair_idx, (criterion, kind) in enumerate(pairs):
        col_slice = slice(2 * pair_idx, 2 * pair_idx + 2)
        role = "primary" if kind == "pass%" and criterion in PRIMARY_CRITERIA else "secondary"

        for is_adjusted, label in enumerate(adjustments):
            est = point[is_adjusted, col_slice]
            sample = draws[:, is_adjusted, col_slice]
            valid = sample[np.isfinite(sample).all(axis=1)]
            delta = valid[:, 0] - valid[:, 1]
            ci = quantile_ci(np.column_stack([valid, delta])) if len(valid) else np.full((2, 3), np.nan)

            rows.append({
                "metric": "GPR" if kind == "pass%" else "mean_gamma",
                "criterion": criterion.replace("_", "/"),
                "role": role,
                "adjustment": label,
                "bootstrap_unit": "patient",
                "n_plans": len(primary),
                "n_patients": primary.patient_id.nunique(),
                "seed": SEED_MAIN,
                "bootstrap_requested": len(draws),
                "bootstrap_valid": len(valid),
                "bootstrap_invalid": len(sample) - len(valid),
                "rho_global": est[0],
                "rho_global_ci_low": ci[0, 0],
                "rho_global_ci_high": ci[1, 0],
                "rho_local": est[1],
                "rho_local_ci_low": ci[0, 1],
                "rho_local_ci_high": ci[1, 1],
                "delta": est[0] - est[1],
                "delta_ci_low": ci[0, 2],
                "delta_ci_high": ci[1, 2],
                "delta_ci_excludes_zero": bool(ci[0, 2] > 0 or ci[1, 2] < 0),
            })
    return rows


def run_associations(raw: pd.DataFrame, primary: pd.DataFrame, output: Path) -> list[dict]:
    """Compute TTDF and secondary metric correlations, paired differences, and CSI sensitivity."""
    patient_clusters = get_patient_clusters(primary.patient_id)
    site_codes = pd.factorize(primary.site_en, sort=False)[0]
    values = primary[["TTDF", *SECONDARY_METRICS, *OUTCOMES]].to_numpy(float)
    primary_indices = [OUTCOMES.index(y) for y in PRIMARY_OUTCOMES]
    primary_values = primary[["TTDF", *PRIMARY_OUTCOMES]].to_numpy(float)
    parameters = primary[PLANNING_PARAMETERS].to_numpy(float)
    site_design = pd.get_dummies(primary.site_en, dtype=float).to_numpy()
    # Exported widths have small numeric variations within each machine setting.
    adjustment_parameters = parameters.copy()
    adjustment_parameters[:, 0] = np.round(adjustment_parameters[:, 0], decimals=2)

    point_base, point_adj = estimate_ranks_and_adjusted(values, site_codes)
    point_planning = np.stack([
        estimate_planning_adjusted(primary_values, site_design, adjustment_parameters, n)
        for n in (3, 4)
    ])

    # All models use the same patient draws and recompute ranks and residuals.
    rng = np.random.default_rng(SEED_MAIN)
    base_draws = np.empty((BOOTSTRAP_DRAWS, len(point_base), len(OUTCOMES)))
    adj_draws = np.empty((BOOTSTRAP_DRAWS, len(OUTCOMES)))
    planning_draws = np.empty((BOOTSTRAP_DRAWS, 2, len(PRIMARY_OUTCOMES)))

    for b in range(BOOTSTRAP_DRAWS):
        idx = resample_clusters(patient_clusters, rng)
        base_draws[b], adj_draws[b] = estimate_ranks_and_adjusted(values[idx], site_codes[idx])
        for j, n in enumerate((3, 4)):
            planning_draws[b, j] = estimate_planning_adjusted(
                primary_values[idx], site_design[idx], adjustment_parameters[idx], n,
            )

    # Save paired global vs local comparison
    combined = np.stack([base_draws[:, 0], adj_draws], axis=1)
    main_rows = build_paired_rows(primary, np.stack([point_base[0], point_adj]), combined,
                                  ADJUSTMENTS[:2], PAIRS)
    save_table(main_rows, output, "ttdf_paired_and_adjusted")

    primary_pairs = [(criterion, kind) for criterion, kind in PAIRS if criterion in PRIMARY_CRITERIA]
    save_table(build_paired_rows(
        primary, np.vstack([point_base[0, primary_indices], point_adj[primary_indices], point_planning]),
        np.concatenate([combined[:, :, primary_indices], planning_draws], axis=1),
        ADJUSTMENTS, primary_pairs,
    ), output, "planning_parameter_associations")
    ranks = rank_and_center(np.column_stack([primary_values[:, 0], parameters]))
    save_table({"parameter": PLANNING_PARAMETERS,
                "rho_with_ttdf": fast_correlation(ranks[:, :1], ranks[:, 1:])[0]},
               output, "planning_parameter_ttdf_correlations")

    # Secondary complexity metrics
    ci_low, ci_high = quantile_ci(base_draws[:, :, primary_indices])
    save_table([
        {
            "metric": metric,
            "outcome": outcome,
            "rho": point_base[i, primary_indices[j]],
            "lo": ci_low[i, j],
            "hi": ci_high[i, j],
            "ci_excludes_zero": bool(ci_low[i, j] > 0 or ci_high[i, j] < 0),
            "bootstrap_requested": BOOTSTRAP_DRAWS,
        }
        for i, metric in enumerate(SECONDARY_METRICS, start=1)
        for j, outcome in enumerate(PRIMARY_OUTCOMES)
    ], output, "secondary_associations")

    # Craniospinal irradiation (CSI) sensitivity analysis
    full_values = raw[["TTDF", *PRIMARY_OUTCOMES]].to_numpy(float)
    full_clusters = get_patient_clusters(raw.patient_id)
    full_point = rank_and_center(full_values)

    rng_site = np.random.default_rng(SEED_SITE)
    full_draws = []
    for _ in range(BOOTSTRAP_DRAWS):
        sampled = rank_and_center(full_values[resample_clusters(full_clusters, rng_site)])
        full_draws.append(fast_correlation(sampled[:, :1], sampled[:, 1:])[0])
    f_lo, f_hi = quantile_ci(np.asarray(full_draws))
    save_table({
        "outcome": PRIMARY_OUTCOMES,
        "primary_rho": point_base[0, primary_indices],
        "primary_lo": ci_low[0],
        "primary_hi": ci_high[0],
        "full_rho": fast_correlation(full_point[:, :1], full_point[:, 1:])[0],
        "full_lo": f_lo,
        "full_hi": f_hi,
        "bootstrap_requested": BOOTSTRAP_DRAWS,
    }, output, "csi_sensitivity")
    return main_rows


# ==============================================================================
# 4. Within-Site Correlations & Diagnostics
# ==============================================================================
def run_site_associations(primary: pd.DataFrame, output: Path) -> None:
    """Compute within-site Spearman rho and leave-one-out sensitivity diagnostics."""
    rows = []

    for site in SITES:
        site_df = primary.loc[primary.site_en == site]
        clusters = get_patient_clusters(site_df.patient_id)
        vals = site_df[["TTDF", GLOBAL_GPR, LOCAL_GPR]].to_numpy(float)

        ranked = rank_and_center(vals)
        point_corr = fast_correlation(ranked[:, :1], ranked[:, 1:])[0]

        # Restart the study seed for each site and keep only valid Global–Local pairs.
        rng = np.random.default_rng(SEED_SITE)
        draws = []
        for _ in range(BOOTSTRAP_DRAWS):
            sample = rank_and_center(vals[resample_clusters(clusters, rng)])
            draws.append(fast_correlation(sample[:, :1], sample[:, 1:])[0])
        draws = np.asarray(draws)
        valid_draws = draws[np.isfinite(draws).all(axis=1)]
        lo, hi = quantile_ci(valid_draws) if len(valid_draws) else np.full((2, 2), np.nan)

        # Leave-one-out sensitivity analysis for global 3%/2mm
        loo = np.asarray([
            stats.spearmanr(*np.delete(vals[:, :2], i, axis=0).T).statistic
            for i in range(len(site_df))
        ])

        for j, outcome in enumerate([GLOBAL_GPR, LOCAL_GPR]):
            rows.append({
                "site": site,
                "outcome": outcome,
                "n": len(site_df),
                "patients": len(clusters),
                "rho": point_corr[j],
                "lo": lo[j],
                "hi": hi[j],
                "bootstrap_requested": BOOTSTRAP_DRAWS,
                "valid_draws": len(valid_draws),
                "loo_low": loo.min() if j == 0 else None,
                "loo_high": loo.max() if j == 0 else None,
                "max_excursion": np.abs(loo - point_corr[0]).max() if j == 0 else None,
            })
    save_table(rows, output, "site_associations")


# ==============================================================================
# 5. Statistical Process Control Limits
# ==============================================================================
def calculate_spc_limits(train: np.ndarray, test: np.ndarray, site: str, outcome: str) -> dict:
    """Calculate Shewhart and weighted-standard-deviation limits."""
    # WSD d2 constants tabulated for Px from 0.50 to 0.70
    p_grid = np.linspace(0.50, 0.70, 11)
    d2_grid = [1.128, 1.126, 1.121, 1.112, 1.099, 1.083, 1.063, 1.039, 1.012, 0.982, 0.947]

    mean = train.mean()
    mr_bar = np.abs(np.diff(train)).mean()
    px = np.mean(train <= mean)
    p_lookup = max(px, 1.0 - px)
    d2 = float(np.interp(p_lookup, p_grid, d2_grid))

    # Classical Shewhart control limits (3-sigma = 2.66 * MRbar)
    shewhart_lcl = mean - 2.66 * mr_bar
    shewhart_ucl = mean + 2.66 * mr_bar

    # Asymmetric WSD control limits
    wsd_lcl = mean - 6.0 * mr_bar * (1.0 - px) / d2
    wsd_ucl = mean + 6.0 * mr_bar * px / d2

    is_gpr = outcome.endswith("pass%")
    return {
        "outcome": outcome,
        "site": site,
        "n": len(train),
        "mean": mean,
        "MRbar": mr_bar,
        "P_x": px,
        "P_lookup": min(p_lookup, 0.70),
        "d2_WSD": d2,
        "d2_endpoint_capped": p_lookup > 0.70,
        "Shewhart_LCL": shewhart_lcl,
        "Shewhart_UCL": shewhart_ucl,
        "WSD_LCL": wsd_lcl,
        "WSD_UCL": wsd_ucl,
        "unfavorable_tail": "lower" if is_gpr else "upper",
        "Shewhart_control_limit": shewhart_lcl if is_gpr else shewhart_ucl,
        "WSD_control_limit": wsd_lcl if is_gpr else wsd_ucl,
        "validation_n": len(test),
        "shewhart_below": int(np.sum(test < shewhart_lcl)),
        "wsd_below": int(np.sum(test < wsd_lcl)),
        "shewhart_above": int(np.sum(test > shewhart_ucl)),
        "wsd_above": int(np.sum(test > wsd_ucl)),
    }


def run_spc_limits(derivation: pd.DataFrame, validation: pd.DataFrame, output: Path) -> None:
    """Compute site-specific Shewhart and WSD limits across all 16 gamma outcomes."""
    records = []
    flags = []
    outcomes = [*PRIMARY_OUTCOMES, *(col for col in OUTCOMES if col not in PRIMARY_OUTCOMES)]
    for outcome in outcomes:
        compare = np.less if outcome.endswith("pass%") else np.greater
        for site in SITES:
            train_vals = derivation.loc[derivation.site_en == site, outcome].to_numpy()
            test = validation.loc[validation.site_en == site]
            test_vals = test[outcome].to_numpy()
            limits = calculate_spc_limits(train_vals, test_vals, site, outcome)
            shewhart = compare(test_vals, limits["Shewhart_control_limit"])
            wsd = compare(test_vals, limits["WSD_control_limit"])
            both = shewhart & wsd
            limits["both_flagged_n"] = int(both.sum())
            records.append(limits)
            flags.append(test[["qa_excel_row", "plan_day", "site_en"]].assign(
                outcome=outcome,
                value=test_vals,
                shewhart_flag=shewhart,
                wsd_flag=wsd,
                both_flag=both,
            ))

    all_df = pd.DataFrame(records)
    save_table(all_df, output, "spc_limits_all_outcomes")

    # Primary global 3%/2mm control limits
    primary_limits = all_df.loc[all_df.outcome == GLOBAL_GPR].drop(columns=["outcome", "both_flagged_n"])
    primary_limits["difference_pp"] = primary_limits.Shewhart_LCL - primary_limits.WSD_LCL
    save_table(primary_limits, output, "spc_limits")

    validation_flags = pd.concat(flags, ignore_index=True)
    save_table(validation_flags, output, "spc_validation_flags")
    reviewed = validation_flags.loc[
        validation_flags.outcome.eq(GLOBAL_GPR)
        & (validation_flags.shewhart_flag | validation_flags.wsd_flag)
    ].merge(validation[["qa_excel_row", "TTDF", "TargetLength", "FieldWidth", "Pitch", "MF"]],
            on="qa_excel_row", validate="one_to_one")
    limit_fields = ["site", "Shewhart_LCL", "WSD_LCL"]
    reviewed = reviewed.merge(primary_limits[limit_fields], left_on="site_en", right_on="site",
                              validate="many_to_one").drop(columns="site")
    reviewed["flag_class"] = np.select(
        [reviewed.both_flag, reviewed.shewhart_flag],
        ["Both", "Shewhart only"],
        default="WSD only",
    )
    save_table(reviewed.sort_values(["both_flag", "plan_day", "qa_excel_row"]),
               output, "spc_global_3pct_2mm_review")
    summary = validation_flags.groupby("outcome", sort=False).agg(
        validation_n=("value", "size"),
        shewhart_flagged_n=("shewhart_flag", "sum"),
        wsd_flagged_n=("wsd_flag", "sum"),
        both_flagged_n=("both_flag", "sum"),
    ).reset_index()
    summary["shewhart_flagged_pct"] = 100 * summary.shewhart_flagged_n / summary.validation_n
    summary["wsd_flagged_pct"] = 100 * summary.wsd_flagged_n / summary.validation_n
    save_table(summary, output, "spc_validation_summary")

    # Overall means and expected means at the observed validation site mix.
    validation_counts = validation.site_en.value_counts()
    derivation_site_means = derivation.groupby("site_en")[OUTCOMES].mean()
    standardized_mean = derivation_site_means.mul(validation_counts, axis=0).sum() / len(validation)
    save_table(pd.DataFrame({
        "derivation_mean": derivation[OUTCOMES].mean(),
        "validation_mean": validation[OUTCOMES].mean(),
        "derivation_mean_at_validation_site_mix": standardized_mean,
    }), output, "spc_temporal_means", index_label="outcome")
    save_table([{
        "subset": name,
        "n_plans": len(subset),
        "first_plan_day": subset.plan_day.min(),
        "last_plan_day": subset.plan_day.max(),
        "plans_on_first_day": int(subset.plan_day.eq(subset.plan_day.min()).sum()),
        "plans_on_last_day": int(subset.plan_day.eq(subset.plan_day.max()).sum()),
    } for name, subset in [("derivation", derivation), ("validation", validation)]],
        output, "chronological_split")


# ==============================================================================
# 6. Chronological Validation
# ==============================================================================
def run_temporal_associations(derivation: pd.DataFrame, validation: pd.DataFrame,
                              main_rows: list[dict], output: Path) -> None:
    """Evaluate temporal stability across 409 derivation and 103 validation plans."""
    subsets = [("derivation", derivation), ("validation", validation)]
    rows = [
        {**row, "subset": "full", "n": row["n_plans"], "patients": row["n_patients"],
         "valid_draws": row["bootstrap_valid"], "primary_spec": True}
        for row in main_rows if row["role"] == "primary" and row["adjustment"] == "unadjusted"
    ]
    for label, subset_df in subsets:
        # Preserve the historical worksheet order for the bootstrap random stream.
        subset_df = subset_df.sort_index().sort_values("plan_day", kind="stable")
        vals = subset_df[["TTDF", *OUTCOMES[:4]]].to_numpy(float)
        ranked = rank_and_center(vals)
        point = fast_correlation(ranked[:, :1], ranked[:, 1:])[0].reshape(2, 2)
        clusters = get_patient_clusters(subset_df.patient_id)

        plan_rng = np.random.default_rng(SEED_TEMPORAL)
        patient_rng = np.random.default_rng(SEED_TEMPORAL)

        samples = np.empty((2, BOOTSTRAP_DRAWS, 2, 2))
        for b in range(BOOTSTRAP_DRAWS):
            # 0: Plan-level bootstrap, 1: Patient-cluster bootstrap
            plan_idx = plan_rng.integers(len(subset_df), size=len(subset_df))
            patient_idx = resample_clusters(clusters, patient_rng)
            for k, idx in enumerate([plan_idx, patient_idx]):
                sampled_ranks = rank_and_center(vals[idx])
                samples[k, b] = fast_correlation(sampled_ranks[:, :1], sampled_ranks[:, 1:])[0].reshape(2, 2)

        for j, criterion in enumerate(PRIMARY_CRITERIA):
            for k, unit in enumerate(["plan", "patient"]):
                draw = samples[k, :, j]
                valid = draw[np.isfinite(draw).all(axis=1)]
                ci = quantile_ci(np.column_stack([valid, valid[:, 0] - valid[:, 1]])) if len(valid) else np.full((2, 3), np.nan)
                rows.append({
                    "subset": label,
                    "n": len(subset_df),
                    "patients": len(clusters),
                    "criterion": criterion.replace("_", "/"),
                    "bootstrap_unit": unit,
                    "seed": SEED_TEMPORAL,
                    "bootstrap_requested": BOOTSTRAP_DRAWS,
                    "valid_draws": len(valid),
                    "primary_spec": unit == "plan",
                    "rho_global": point[j, 0],
                    "rho_global_ci_low": ci[0, 0],
                    "rho_global_ci_high": ci[1, 0],
                    "rho_local": point[j, 1],
                    "rho_local_ci_low": ci[0, 1],
                    "rho_local_ci_high": ci[1, 1],
                    "delta": point[j, 0] - point[j, 1],
                    "delta_ci_low": ci[0, 2],
                    "delta_ci_high": ci[1, 2],
                })
    columns = ["subset", "n", "patients", "criterion", "bootstrap_unit", "seed",
               "bootstrap_requested", "valid_draws", "primary_spec", "rho_global", "rho_global_ci_low",
               "rho_global_ci_high", "rho_local", "rho_local_ci_low", "rho_local_ci_high", "delta",
               "delta_ci_low", "delta_ci_high"]
    save_table(pd.DataFrame(rows)[columns], output, "temporal_associations")


# ==============================================================================
# CLI Entrypoint
# ==============================================================================
def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "analysis_cohort.csv")
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    args, unknown = parser.parse_known_args()

    # Handle positional arguments for legacy invocation
    cohort_path = args.data if not unknown else Path(unknown[0])
    output_path = args.out if len(unknown) <= 1 else Path(unknown[1])
    output_path.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(cohort_path, float_precision="round_trip")
    primary = raw.loc[raw.analysis_cohort == "Primary standard-site"]
    derivation, validation = chronological_split(primary)

    # Run analysis stages
    run_clustering(primary, output_path)
    run_distributions(primary, output_path)
    main_rows = run_associations(raw, primary, output_path)
    run_site_associations(primary, output_path)
    run_spc_limits(derivation, validation, output_path)
    run_temporal_associations(derivation, validation, main_rows, output_path)


if __name__ == "__main__":
    main()
