# Helical Tomotherapy PSQA Plan Complexity and SPC Analysis Package

Analysis code and input-data specification for evaluating helical tomotherapy (HT) plan-complexity metrics, global versus local gamma outcomes, and site-specific statistical process control (SPC) limits (associated study currently under peer review).

## Reproduce the analysis

The analysis was verified with CPython 3.11.7 on Windows x64. Neither the institutional cohort CSV nor pre-generated outputs are bundled in this repository. After obtaining or preparing `data/analysis_cohort.csv` (see [data/README.md](data/README.md)), create a virtual environment and run:

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe run.py
```

On macOS/Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run.py
```

`run.py` reads `data/analysis_cohort.csv` (or `--data PATH`) and executes all statistical analyses, SPC and normality audits, logistic discrimination modeling, plan-characteristics summaries, and summary figure generation. Excel workbooks (`.xlsx`) are written to `results/` (or `--out PATH`), and figures are saved as vector PDF and 300-DPI PNG in `results/figures/`.

The same `--data` and `--out` options are supported by each standalone script:
- `analysis/run_analysis.py` — metric clustering, gamma distributions, TTDF/secondary associations, planning-parameter adjustments, SPC limits, and chronological validation
- `analysis/analyze_ad_norm.py` — Anderson–Darling normality diagnostics across all 16 gamma outcomes
- `analysis/discrimination.py` — exploratory logistic discrimination modeling and grouped cross-validation
- `analysis/plan_characteristics.py` — treatment-plan delivery and geometric summaries
- `analysis/figures.py` — summary figures (also accepts `--results PATH`)

## Input data availability and preparation

Clinical cohort data (`data/analysis_cohort.csv`) are not stored in this repository and can be shared upon reasonable request. See [data/README.md](data/README.md) for the complete column schema and units.

To export a de-identified CSV from an institutional study workbook (`Analysis_Cohort` worksheet):

```powershell
.venv\Scripts\python.exe tools/export_cohort.py --workbook PATH_TO_WORKBOOK.xlsx
.venv\Scripts\python.exe run.py
```

The exporter preserves worksheet order, assigns sequential patient pseudonyms (`P0001`, `P0002`, ...), converts plan dates to relative elapsed study days (`plan_day`), and writes checksum metadata to `data/cohort_summary.json` (both `data/*.csv` and `data/*.json` are git-ignored).

## Analysis modules and generated outputs

Running `run.py` populates `results/` with the following topic-organized workbooks and figures:

| Analysis module | Generated files in `results/` |
|---|---|
| Treatment-plan characteristics | `plan_characteristics.xlsx` |
| Complexity-metric correlation and clustering | `spearman_matrix.xlsx`, `clusters.xlsx`, `linkage.xlsx`, `figures/complexity_heatmap.{pdf,png}` |
| Pooled and site-specific gamma distributions | `gamma_distributions.xlsx`, `figures/site_distributions.{pdf,png}` |
| Primary TTDF associations and paired global–local differences | `ttdf_paired_and_adjusted.xlsx`, `figures/ttdf_associations.{pdf,png}`, `figures/criterion_sensitivity.{pdf,png}` |
| Secondary complexity-metric associations | `secondary_associations.xlsx` |
| Within-site TTDF–GPR associations and LOO diagnostics | `site_associations.xlsx`, `figures/site_associations.{pdf,png}` |
| CSI expanded-cohort sensitivity analysis | `csi_sensitivity.xlsx` |
| Planning-parameter partial-rank adjustments | `planning_parameter_associations.xlsx`, `planning_parameter_ttdf_correlations.xlsx` |
| Site-specific Shewhart and WSD control limits & validation flags | `spc_limits.xlsx`, `spc_limits_all_outcomes.xlsx`, `spc_validation_summary.xlsx`, `spc_validation_flags.xlsx`, `spc_global_3pct_2mm_review.xlsx` |
| Chronological split, temporal means, and temporal stability | `chronological_split.xlsx`, `spc_temporal_means.xlsx`, `temporal_associations.xlsx` |
| Anderson–Darling normality diagnostics | `ad_normality.xlsx`, `ad_normality_summary.xlsx` |
| Exploratory logistic discrimination (global 2%/2 mm GPR <95%) | `discrimination_events.xlsx`, `discrimination_models.xlsx`, `discrimination_coefficients.xlsx`, `discrimination_cv_repeats.xlsx`, `discrimination_standardization.xlsx` |

## Analysis settings

- **Metric clustering.** Complete linkage on `1 − |ρ|` (Spearman rank correlation), cut at height 0.50. Retained representative metrics are TTDF, CFNS75, L0NS, MI, and centroid.
- **Associations.** All bootstrap analyses use 3,000 draws and percentile 95% confidence intervals (CIs). Spearman correlations recompute average ranks within each draw. Main, within-site, and expanded-cohort analyses resample patient pseudo-IDs with replacement and retain all plans per sampled ID; global and local outcomes share each draw. Partially adjusted models recompute residuals within every draw. Secondary complexity-metric correlations report unadjusted pointwise CIs across 32 exploratory pairs. Seeds: `20260929` (main analyses), `20260930` (within-site and expanded-cohort analyses), and `20261116` (temporal validation and discrimination bootstrap). Within-site analyses exclude draws with undefined correlations for either normalization and record valid draw counts.
- **Planning-parameter adjustment.** TTDF and gamma ranks are residualized on treatment-site indicators and ranked field width, pitch, and MF, followed by ranked target length. Field width is rounded to two decimal places before ranking to group the three nominal jaw settings; descriptive TTDF–field-width correlations use unrounded exported values.
- **SPC and temporal analyses.** Control limits use Shewhart I-chart (`X̄ ± 2.66 × MR̄`) and Chang–Bai weighted-standard-deviation (WSD) methods, with `d2` linearly interpolated over `P_lookup` in `[0.50, 0.70]` and endpoint-clamped without extrapolation (`d2_endpoint_capped` flags clamped cases). Plans are ordered by `(plan_day, qa_excel_row)` and split into 409 derivation and 103 validation plans without dividing tied days. In `temporal_associations.xlsx`, full-cohort rows reuse the main paired patient-cluster bootstrap; derivation/validation subset rows report paired plan-level bootstrap (`primary_spec = TRUE`, preserving original worksheet order within day) alongside paired patient-cluster bootstrap results for comparison.
- **SPC flags.** Evaluated across all 16 gamma outcomes (4 criteria × 2 normalizations × 2 summary statistics). Unfavorable exceedances are strictly below the LCL for GPR or strictly above the UCL for mean gamma.
- **Normality.** Uses `statsmodels.stats.diagnostic.normal_ad` with estimated mean and variance (`p < 0.05`) alongside SciPy interpolated Anderson–Darling statistics for comparison. Site-test summaries count the 96 site-specific tests and exclude the 16 pooled tests.
- **Discrimination.** Evaluates unpenalized logistic regression (`C = inf`, requiring `scikit-learn >= 1.8` as pinned in `requirements.txt`) for global 2%/2 mm GPR <95%. Apparent AUCs and ΔAUC use 3,000 patient-cluster bootstrap draws of fixed fitted scores. Grouped cross-validation reports 20 repeats of 5-fold `StratifiedGroupKFold` with z-standardization fitted inside each training fold. Chronological validation uses the 409/103 split and records patient overlap.
- **Treatment-plan characteristics.** Dose per fraction is reconstructed from `TreatmentTime / TTDF / 100` (Gy) and nominal prescription as `dose_per_fraction * fractions` across plans with recorded fraction counts.

## Optional: Extract plan complexity metrics from DICOM RTPLAN

`extraction/` contains a local adaptation of **TCoMX v1.0** (Cavinato et al., *J Appl Clin Med Phys* 2023;24(1):e13781; GNU GPLv3) together with the batch wrapper `extraction/extract_metrics.m`. File integrity digests are recorded in `extraction/SHA256SUMS`. Note that `extraction/libs/TCoMX.m` is retained as an upstream reference script; use `extraction/extract_metrics.m` as the entrypoint.

In MATLAB (with Image Processing Toolbox), prepare a CSV manifest with columns `plan_key`, `rtplan_path`, and `clinical_dose_per_fraction_gy` (for multi-prescription SIB plans, provide the highest prescription level's dose per fraction; if absent, TTDF is left `NaN` while all other metrics are still computed):

```matlab
addpath('extraction');
extract_metrics('manifest.csv', 'extracted_metrics.csv');
```

Running the Python analysis pipeline (`run.py`) requires only a formatted cohort CSV (`data/analysis_cohort.csv`) and does not require MATLAB or access to institutional DICOM RTPLAN files.

## License and citation

- **Analysis and export code** (`run.py`, `analysis/`, `tools/`): [MIT License](LICENSE).
- **MATLAB complexity extraction code** (`extraction/`): [GNU GPLv3](extraction/LICENSE), adapted from TCoMX v1.0.

Citation metadata for this software package is provided in [CITATION.cff](CITATION.cff). Publication reference details will be added upon article acceptance.
