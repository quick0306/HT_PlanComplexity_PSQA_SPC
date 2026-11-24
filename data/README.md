# Input Data Specification

This directory defines the expected CSV input format (`analysis_cohort.csv`) for the helical tomotherapy patient-specific quality assurance (PSQA) plan-complexity and statistical process control (SPC) analysis pipeline.

## Data Availability

The study cohort dataset is **not bundled** in this public repository and is available upon reasonable request. To run the analysis locally, place a formatted cohort CSV at `data/analysis_cohort.csv` (or pass `--data PATH` to `run.py`), or generate it from an institutional study workbook using `tools/export_cohort.py`.

## Export and De-identification Workflow

`tools/export_cohort.py` extracts the required columns from an institutional workbook (`Analysis_Cohort` worksheet), replaces patient identifiers with sequential pseudonyms (`P0001`, `P0002`, ...), converts RT plan dates to elapsed study days (`plan_day`, starting at Day 0) while preserving original QA record ordering (`qa_excel_row`), and writes `data/analysis_cohort.csv` together with integrity metadata (`data/cohort_summary.json`).

## Expected Schema (`analysis_cohort.csv`)

| Column | Type | Unit | Description |
|---|---|---|---|
| `qa_excel_row` | Integer | - | Dense rank of the original QA record order; preserves chronological tie-breaking. |
| `patient_id` | String | - | De-identified patient pseudonym (e.g., `P0001`, `P0002`, ...). Preserves patient clustering for bootstrap resampling. |
| `site_en` | String | - | Anatomical treatment site (`Brain`, `Head & neck`, `Thorax`, `Breast`, `Abdomen`, `Pelvis`, or `CSI`). |
| `analysis_cohort` | String | - | Cohort classification: `Primary standard-site` or `CSI sensitivity`. |
| `plan_day` | Integer | Days | Elapsed days from the earliest treatment plan in the cohort (Day 0). |
| `TTDF` | Float | s/cGy | Treatment time per dose fraction, normalized to fraction dose. |
| `MF` | Float | - | Modulation factor: ratio of maximum to mean non-zero leaf open times. |
| `mFLOT` | Float | - | Mean fractional leaf-open time. |
| `CFNS75` | Float | - | Cumulative fraction of positive fractional leaf-open times below 0.75. |
| `L0NS` | Float | - | Fraction of open leaves with no open neighbors. |
| `CLS` | Float | - | Closed-leaf score. |
| `fDISC` | Float | - | Fraction of discontinuous projections in the sinogram. |
| `centroid` | Float | - | Longitudinal centroid of sinogram openings. |
| `PSTV` | Float | - | Plan sinogram time variation. |
| `MI` | Float | - | Modulation index. |
| `TreatmentTime` | Float | s | Planned delivery time in seconds extracted from DICOM RTPLAN. |
| `TargetLength` | Float | mm | Longitudinal target length (couch translation minus computed field width). |
| `FieldWidth` | Float | mm | Effective jaw field width derived from RTPLAN delivery geometry (approximately 10.48, 25.12, or 50.48 mm, corresponding to nominal 10, 25, and 50 mm jaw settings). |
| `Pitch` | Float | - | Nominal helical pitch. |
| `fractions` | Float | - | Prescribed fraction count (may be blank when total fraction count is unavailable). |
| `global_*_pass%` | Float | % | Global gamma passing rate (0–100%) for criteria `3%_3mm`, `3%_2mm`, `2%_2mm`, `2%_1mm`. |
| `local_*_pass%` | Float | % | Local gamma passing rate (0–100%) for criteria `3%_3mm`, `3%_2mm`, `2%_2mm`, `2%_1mm`. |
| `global_*_mean_gamma` | Float | - | Global mean gamma index for criteria `3%_3mm`, `3%_2mm`, `2%_2mm`, `2%_1mm`. |
| `local_*_mean_gamma` | Float | - | Local mean gamma index for criteria `3%_3mm`, `3%_2mm`, `2%_2mm`, `2%_1mm`. |
