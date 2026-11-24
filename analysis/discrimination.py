"""Exploratory logistic discrimination analysis.

Models the event global 2%/2 mm GPR < 95%,
using clinical site and complexity metrics.

Outputs:
- discrimination_models.xlsx, discrimination_coefficients.xlsx
- Site-stratified event rates (discrimination_events.xlsx)
- 20-repeat 5-fold grouped cross-validation AUC distributions (discrimination_cv_repeats.xlsx)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
import statsmodels.api as sm

if __package__:
    from .result_tables import save_table
    from .study_settings import BOOTSTRAP_DRAWS, chronological_split
else:
    from result_tables import save_table
    from study_settings import BOOTSTRAP_DRAWS, chronological_split


ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = "global_2%_2mm_pass%"
SITES = ["Abdomen", "Brain", "Breast", "Head & neck", "Pelvis", "Thorax"]
METRICS = ["TTDF", "CFNS75", "L0NS", "MI", "centroid"]

# Nested multivariable logistic regression model specifications
MODELS = {
    "Site only": [],
    "Site + TTDF": ["TTDF"],
    "Site + TTDF + L0NS + centroid": ["TTDF", "L0NS", "centroid"],
    "Site + all five metrics": METRICS,
}


def fit_and_predict_risk(site: np.ndarray, metrics: pd.DataFrame, event: np.ndarray,
                         train: np.ndarray, test: np.ndarray) -> dict[str, np.ndarray]:
    """Standardize metrics on train set and compute predicted event probabilities on test set."""
    train_metrics = metrics.iloc[train]
    train_sd = train_metrics.std(ddof=1).replace(0.0, 1.0)
    scaled = (metrics - train_metrics.mean()) / train_sd
    predictions = {}
    for model_name, terms in MODELS.items():
        x_train = np.column_stack([site[train], scaled.iloc[train][terms]])
        x_test = np.column_stack([site[test], scaled.iloc[test][terms]])
        clf = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=2000, tol=1e-8)
        clf.fit(x_train, event[train])
        predictions[model_name] = clf.predict_proba(x_test)[:, 1]
    return predictions


def bootstrap_auc_ci(event: np.ndarray, scores: np.ndarray,
                      patient_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Estimate 95% bootstrap CIs for AUC and delta-AUC relative to Site-only baseline.

    Uses 3,000 patient-cluster resamples and Mann-Whitney U ranking for computational speed.
    """
    clusters = [np.flatnonzero(patient_ids == pid) for pid in np.unique(patient_ids)]
    rng = np.random.default_rng(20261116)
    draws = []

    for _ in range(BOOTSTRAP_DRAWS):
        selected_clusters = rng.integers(len(clusters), size=len(clusters))
        rows = np.concatenate([clusters[i] for i in selected_clusters])
        ev = event[rows].astype(bool)
        n_pos = ev.sum()
        n_neg = len(ev) - n_pos
        if n_pos > 0 and n_neg > 0:
            # Fast AUC calculation via rank sums
            ranks = rankdata(scores[rows], axis=0)
            u_stat = ranks[ev].sum(axis=0) - n_pos * (n_pos + 1) / 2.0
            draws.append(u_stat / (n_pos * n_neg))

    samples = np.asarray(draws)
    auc_ci = np.quantile(samples, [0.025, 0.975], axis=0).T
    delta_ci = np.quantile(samples - samples[:, :1], [0.025, 0.975], axis=0).T
    return auc_ci, delta_ci, len(samples)


def run_grouped_cv(site: np.ndarray, metrics: pd.DataFrame, event: np.ndarray,
                   patient_ids: np.ndarray) -> pd.DataFrame:
    """Run 20 repeats of 5-fold Stratified Group K-Fold cross-validation."""
    records = []
    for repeat in range(20):
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=1102 + repeat)
        oof_preds = {m: np.empty(len(event)) for m in MODELS}

        for train_idx, test_idx in cv.split(site, event, patient_ids):
            preds = fit_and_predict_risk(site, metrics, event, train_idx, test_idx)
            for m, p in preds.items():
                oof_preds[m][test_idx] = p

        records.append({m: roc_auc_score(event, pred) for m, pred in oof_preds.items()})
    return pd.DataFrame(records, index=pd.RangeIndex(1, 21, name="repeat"))


def discriminate(data_path: Path, output_dir: Path) -> None:
    """Fit full cohort models, compute statistics, cross-validation, and temporal validation."""
    cohort = (
        pd.read_csv(data_path, float_precision="round_trip")
        .loc[lambda df: df.analysis_cohort.eq("Primary standard-site")]
        .reset_index(drop=True)
    )

    # Event definition: Global 2%/2mm passing rate < 95%
    event = cohort[ENDPOINT].lt(95).to_numpy(dtype=int)
    patient = cohort.patient_id.to_numpy(dtype=str)
    metrics = cohort[METRICS]
    metric_mean = metrics.mean()
    metric_sd = metrics.std(ddof=1)
    standardized = (metrics - metric_mean) / metric_sd

    # Site dummy variables (Abdomen as reference site)
    site_dummies = np.column_stack([cohort.site_en.eq(name) for name in SITES[1:]]).astype(float)
    site_labels = ["Intercept", *(f"{s} vs Abdomen" for s in SITES[1:])]

    # Site event summary
    event_counts = (
        cohort.assign(events=event)
        .groupby("site_en")
        .events.agg(n="size", events="sum")
        .reindex(SITES)
    )
    event_counts.loc["Overall"] = event_counts.sum()
    event_counts["nonevents"] = event_counts.n - event_counts.events
    event_counts["event_pct"] = event_counts.events * 100.0 / event_counts.n

    # Fit GLM models with patient-clustered robust standard errors
    model_summaries, coefficient_dfs, scores = [], [], {}
    for model_name, terms in MODELS.items():
        design = np.column_stack([np.ones(len(event)), site_dummies, standardized[terms]])
        fit = sm.GLM(event, design, family=sm.families.Binomial()).fit(
            maxiter=200, cov_type="cluster", cov_kwds={"groups": patient}
        )
        scores[model_name] = fit.predict()
        model_summaries.append({
            "model": model_name,
            "log_likelihood": fit.llf,
            "parameters": design.shape[1],
            "aic": fit.aic,
            "deviance": fit.deviance,
            "converged": fit.converged,
        })
        ci_lo, ci_hi = fit.conf_int().T
        coefficient_dfs.append(pd.DataFrame({
            "model": model_name,
            "predictor": site_labels + terms,
            "beta": fit.params,
            "cluster_se": fit.bse,
            "wald_p": fit.pvalues,
            "beta_lo": ci_lo,
            "beta_hi": ci_hi,
            "odds_ratio": np.exp(fit.params),
            "or_lo": np.exp(ci_lo),
            "or_hi": np.exp(ci_hi),
        }))

    # Univariable metric models
    for metric in ("TTDF", "L0NS"):
        m_name = f"{metric} univariable"
        scores[m_name] = metrics[metric].to_numpy()
        model_summaries.append({"model": m_name})

    # Likelihood ratio tests relative to Site-only baseline
    res_df = pd.DataFrame(model_summaries).set_index("model")
    baseline_ll = res_df.loc["Site only", "log_likelihood"]
    baseline_df = res_df.loc["Site only", "parameters"]
    res_df["lr_stat_vs_site"] = 2 * (res_df.log_likelihood - baseline_ll)
    res_df["lr_df_vs_site"] = res_df.parameters - baseline_df
    res_df["lr_p_vs_site"] = chi2.sf(res_df.lr_stat_vs_site, res_df.lr_df_vs_site)
    res_df.loc["Site only", ["lr_stat_vs_site", "lr_df_vs_site", "lr_p_vs_site"]] = np.nan

    # Compare the full model with the site-plus-TTDF model.
    full = res_df.loc["Site + all five metrics"]
    site_ttdf = res_df.loc["Site + TTDF"]
    lr = 2 * (full.log_likelihood - site_ttdf.log_likelihood)
    degrees = full.parameters - site_ttdf.parameters
    res_df.loc["Site + all five metrics", ["lr_stat_vs_site_ttdf", "lr_df_vs_site_ttdf", "lr_p_vs_site_ttdf"]] = (
        lr, degrees, chi2.sf(lr, degrees),
    )

    # Apparent AUC and bootstrap confidence intervals
    res_df["apparent_auc"] = [roc_auc_score(event, sc) for sc in scores.values()]
    auc_ci, delta_ci, valid_draws = bootstrap_auc_ci(event, np.column_stack(list(scores.values())), patient)
    res_df[["auc_lo", "auc_hi"]] = auc_ci
    res_df["n_plans"] = len(cohort)
    res_df["n_patients"] = cohort.patient_id.nunique()
    res_df["n_events"] = int(event.sum())
    res_df["auc_bootstrap_requested"] = BOOTSTRAP_DRAWS
    res_df["auc_bootstrap_valid"] = valid_draws
    res_df["delta_auc_vs_site"] = res_df.apparent_auc - res_df.apparent_auc.iloc[0]
    res_df[["delta_lo", "delta_hi"]] = delta_ci
    res_df.loc[["Site only", "TTDF univariable", "L0NS univariable"],
               ["delta_auc_vs_site", "delta_lo", "delta_hi"]] = np.nan

    # 20-repeat 5-fold cross-validation
    cv_df = run_grouped_cv(site_dummies, metrics, event, patient)
    res_df["cv_median"] = cv_df.median()
    res_df["cv_min"] = cv_df.min()
    res_df["cv_max"] = cv_df.max()
    res_df["cv_q1"] = cv_df.quantile(0.25)
    res_df["cv_q3"] = cv_df.quantile(0.75)
    cv_delta = cv_df.drop(columns="Site only").sub(cv_df["Site only"], axis=0)
    res_df["cv_delta_median"] = cv_delta.median()
    res_df["cv_delta_min"] = cv_delta.min()
    res_df["cv_delta_max"] = cv_delta.max()

    # Chronological validation (train on first 409 derivation plans, test on last 103)
    derivation, validation = chronological_split(cohort)
    train_idx, test_idx = derivation.index.to_numpy(), validation.index.to_numpy()
    temporal_preds = fit_and_predict_risk(site_dummies, metrics, event, train_idx, test_idx)
    res_df["temporal_auc"] = pd.Series({m: roc_auc_score(event[test_idx], p) for m, p in temporal_preds.items()})
    res_df.loc[list(MODELS), "temporal_n"] = len(test_idx)
    res_df.loc[list(MODELS), "temporal_events"] = event[test_idx].sum()
    res_df.loc[list(MODELS), "temporal_overlap_ids"] = len(set(patient[train_idx]) & set(patient[test_idx]))

    # Save output tables
    output_dir.mkdir(parents=True, exist_ok=True)
    save_table(event_counts, output_dir, "discrimination_events", index_label="site")
    save_table(res_df, output_dir, "discrimination_models", index_label="model")
    save_table(pd.concat(coefficient_dfs, ignore_index=True), output_dir, "discrimination_coefficients")
    cv_audit = cv_df.join(cv_delta.add_suffix("_delta_vs_site"))
    save_table(cv_audit, output_dir, "discrimination_cv_repeats", index_label="repeat_index")
    save_table(pd.DataFrame({"mean": metric_mean, "sample_sd": metric_sd}),
               output_dir, "discrimination_standardization", index_label="metric")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "analysis_cohort.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    discriminate(args.data, args.out)


if __name__ == "__main__":
    main()
