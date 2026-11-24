"""Shared settings and chronological split helper for HT PSQA analyses."""

import pandas as pd


BOOTSTRAP_DRAWS = 3000
DERIVATION_PLANS = 409


def chronological_split(primary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Order by source plan day and QA record sequence, keeping cutoff-day plans together."""
    if primary.plan_day.isna().any() or primary.qa_excel_row.isna().any():
        raise ValueError("Chronological analysis requires a plan day and QA sequence for every plan.")
    ordered = primary.sort_values(["plan_day", "qa_excel_row"])
    derivation = ordered.iloc[:DERIVATION_PLANS]
    validation = ordered.iloc[DERIVATION_PLANS:]
    if validation.empty or derivation.plan_day.iloc[-1] >= validation.plan_day.iloc[0]:
        raise ValueError("The 409-plan split requires a later, nonempty validation cohort and no split dates.")
    return derivation, validation
