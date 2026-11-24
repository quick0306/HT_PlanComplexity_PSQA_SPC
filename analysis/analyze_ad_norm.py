"""Anderson–Darling normality diagnostics for all 16 gamma outcomes."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.diagnostic import normal_ad

if __package__:
    from .result_tables import save_table
    from .run_analysis import OUTCOMES, PRIMARY_CRITERIA, SITES
else:
    from result_tables import save_table
    from run_analysis import OUTCOMES, PRIMARY_CRITERIA, SITES


ROOT = Path(__file__).resolve().parents[1]
ALPHA = 0.05


def run_normality(primary: pd.DataFrame, output: Path) -> None:
    """Write pooled and site-specific diagnostics and counts of site tests."""
    groups = [("Overall", primary), *((site, primary.loc[primary.site_en == site]) for site in SITES)]
    rows = []
    for outcome in OUTCOMES:
        mode, dose, distance, measure = outcome.split("_", maxsplit=3)
        criterion = f"{dose}/{distance}"
        for site, group in groups:
            values = group[outcome].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            unique = np.unique(values).size
            ad, p_value = normal_ad(values) if unique > 1 else (np.nan, np.nan)
            check = stats.anderson(values, dist="norm", method="interpolate") if unique > 1 else None
            rows.append({
                "criterion": criterion,
                "mode": mode,
                "outcome": "GPR" if measure == "pass%" else "Mean gamma",
                "site": site,
                "n": len(values),
                "AD_A2": ad,
                "p_normal_ad": p_value,
                "reject_0_05": bool(p_value < ALPHA) if np.isfinite(p_value) else None,
                "scipy_AD_A2": check.statistic if check is not None else np.nan,
                "scipy_p_interpolate": check.pvalue if check is not None else np.nan,
                "skewness": stats.skew(values, bias=False) if unique > 1 else np.nan,
                "min": values.min() if len(values) else np.nan,
                "max": values.max() if len(values) else np.nan,
                "number_unique": unique,
                "sd": values.std(ddof=1) if len(values) > 1 else np.nan,
            })
    detail = pd.DataFrame(rows)
    save_table(detail, output, "ad_normality")

    site_tests = detail.loc[detail.site != "Overall"]
    primary_criteria = [c.replace("_", "/") for c in PRIMARY_CRITERIA]
    scopes = {
        **dict(tuple(site_tests.groupby("criterion", sort=False))),
        "Primary criteria": site_tests.loc[site_tests.criterion.isin(primary_criteria)],
        "All criteria": site_tests,
    }
    summary = []
    for scope, tests in scopes.items():
        valid = tests.reject_0_05.dropna()
        summary.append({
            "site_test_scope": scope,
            "valid_tests": len(valid),
            "rejected_0_05": int(valid.sum()),
            "not_rejected_0_05": int(len(valid) - valid.sum()),
            "undefined_tests": len(tests) - len(valid),
        })
    save_table(summary, output, "ad_normality_summary")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "analysis_cohort.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()

    cohort = pd.read_csv(args.data, float_precision="round_trip")
    primary = cohort.loc[cohort.analysis_cohort.eq("Primary standard-site")]
    args.out.mkdir(parents=True, exist_ok=True)
    run_normality(primary, args.out)


if __name__ == "__main__":
    main()
