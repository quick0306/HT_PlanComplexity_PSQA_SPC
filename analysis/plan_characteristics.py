"""Summarize treatment plan delivery and geometric characteristics.

Outputs:
- plan_characteristics.xlsx: prescription, fractionation, target length, field width,
  pitch, TTDF, and beam-on time across six anatomical sites and the pooled cohort.
"""

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

if __package__:
    from .result_tables import save_table
else:
    from result_tables import save_table


ROOT = Path(__file__).resolve().parents[1]
METRICS = ["prescription_gy", "dose_per_fraction_gy", "fractions", "TargetLength",
           "FieldWidth", "Pitch", "TTDF", "TreatmentTime"]
SITE_ORDER = ["Brain", "Head & neck", "Thorax", "Breast", "Abdomen", "Pelvis"]


def summarize_plan_characteristics(cohort: pd.DataFrame) -> pd.DataFrame:
    """Calculate median and range for delivery parameters stratified by treatment site."""
    # The curated workbook stores delivery time, TTDF, and audited fraction counts.
    dose_per_fraction = cohort.TreatmentTime / cohort.TTDF.where(cohort.TTDF > 0) / 100
    cohort = cohort.assign(dose_per_fraction_gy=dose_per_fraction,
                           prescription_gy=dose_per_fraction * cohort.fractions)
    groups = [("Overall", cohort)] + [(site, cohort[cohort.site_en == site]) for site in SITE_ORDER]
    rows = []
    for site, group in groups:
        for metric in METRICS:
            values = group[metric]
            q1, median, q3 = values.quantile([0.25, 0.5, 0.75])
            rows.append({
                "site": site,
                "metric": metric,
                "n_plans": len(group),
                "n_available": values.count(),
                "median": median,
                "q1": None if metric == "fractions" else q1,
                "q3": None if metric == "fractions" else q3,
                "minimum": values.min(),
                "maximum": values.max(),
            })
    return pd.DataFrame(rows)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "analysis_cohort.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()

    cohort = pd.read_csv(args.data, float_precision="round_trip")
    primary = cohort[cohort.analysis_cohort == "Primary standard-site"]

    args.out.mkdir(parents=True, exist_ok=True)
    summary_df = summarize_plan_characteristics(primary)
    save_table(summary_df, args.out, "plan_characteristics")


if __name__ == "__main__":
    main()
