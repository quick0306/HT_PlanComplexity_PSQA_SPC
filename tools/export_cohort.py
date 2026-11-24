"""Export the analysis columns from the study workbook in Data to the BPEX package.

Reads the frozen study workbook from 'Data', anonymizes patient identifiers,
normalizes RT plan dates to relative study days, and exports the analysis CSV.
"""

from argparse import ArgumentParser
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path

import pandas as pd


# Directory paths
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKBOOK = PACKAGE_ROOT.parent / "Data" / "TOMO_Gamma_analysis_simplified_with_patient_info_20260426.xlsx"
DEFAULT_OUT = PACKAGE_ROOT / "data"

# Ten plan complexity metrics from TCoMX extraction
METRICS = [
    "TTDF", "MF", "mFLOT", "CFNS75", "L0NS",
    "CLS", "fDISC", "centroid", "PSTV", "MI"
]

# Plan delivery and geometric parameters
PLAN_FIELDS = ["TreatmentTime", "TargetLength", "FieldWidth", "Pitch", "fractions"]

# Gamma passing rate and mean gamma columns across 4 criteria (global & local)
GAMMA = [
    f"{mode}_{criterion}_{outcome}"
    for mode in ["global", "local"]
    for criterion in ["3%_3mm", "3%_2mm", "2%_2mm", "2%_1mm"]
    for outcome in ["pass%", "mean_gamma"]
]

SHEET = "Analysis_Cohort"


def export(workbook: Path, output: Path) -> None:
    """Extract and anonymize analysis cohort from source Excel workbook."""
    source_fields = ["qa_excel_row", "patient_id", "site_en", "analysis_cohort", "RTPlanDate_iso"]
    numeric_fields = METRICS + PLAN_FIELDS + GAMMA
    workbook_bytes = workbook.read_bytes()
    source = pd.read_excel(BytesIO(workbook_bytes), sheet_name=SHEET, usecols=source_fields + numeric_fields)
    if source.patient_id.isna().any():
        raise ValueError("Every plan needs a source patient ID before de-identification.")

    # Anonymize patient IDs to sequential pseudo-identifiers (P0001, P0002, ...)
    id_codes = pd.factorize(source.patient_id.astype(str), sort=True)[0] + 1
    dates = pd.to_datetime(source.RTPlanDate_iso).dt.normalize()

    # Build de-identified cohort dataframe preserving row ordering and within-day sequence
    cohort = pd.DataFrame({
        "qa_excel_row": source.qa_excel_row.rank(method="dense").astype(int),
        "patient_id": [f"P{code:04d}" for code in id_codes],
        "site_en": source.site_en,
        "analysis_cohort": source.analysis_cohort,
        "plan_day": (dates - dates.min()).dt.days,
        **{col: source[col] for col in numeric_fields},
    })

    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "analysis_cohort.csv"
    cohort.to_csv(csv_path, index=False, lineterminator="\n")

    # Generate summary metadata
    primary = cohort[cohort.analysis_cohort == "Primary standard-site"]
    csi = cohort[cohort.analysis_cohort == "CSI sensitivity"]
    summary = {
        "source_workbook": workbook.name,
        "source_sheet": SHEET,
        "source_sha256": sha256(workbook_bytes).hexdigest(),
        "analysis_cohort_sha256": sha256(csv_path.read_bytes()).hexdigest(),
        "counts": {
            label: {"plans": len(data), "patient_ids": int(data.patient_id.nunique())}
            for label, data in [("primary", primary), ("csi", csi), ("expanded", cohort)]
        },
        "patient_ids_in_both_cohorts": len(set(primary.patient_id) & set(csi.patient_id)),
        "primary_site_counts": primary.site_en.value_counts().sort_index().to_dict(),
        "missing_fractions": {
            "primary": int(primary.fractions.isna().sum()),
            "expanded": int(cohort.fractions.isna().sum()),
        },
    }
    (output / "cohort_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Exported {len(cohort)} plans ({len(primary)} primary, {len(csi)} CSI) to {csv_path}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK, help="Path to raw Excel workbook in Data/")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output directory for CSV and metadata")
    args = parser.parse_args()
    export(args.workbook, args.out)
