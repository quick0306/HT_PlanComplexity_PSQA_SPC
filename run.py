"""Reproduce HT PSQA statistical analyses, audits, and summary figures."""

import argparse
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "analysis_cohort.csv",
                        help="Path to curated analysis cohort CSV")
    parser.add_argument("--out", type=Path, default=ROOT / "results",
                        help="Target directory for generated tables and figures")
    args = parser.parse_args()

    data_path = args.data.resolve()
    if not data_path.is_file():
        parser.error(
            f"Input cohort CSV not found: {data_path}\n"
            "Place the analysis cohort CSV at data/analysis_cohort.csv (see data/README.md for schema) "
            "or specify another location with --data PATH."
        )
    out_dir = args.out.resolve()
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    pipeline = [
        ("Primary statistical analysis & SPC audit", "run_analysis.py",
         ["--data", str(data_path), "--out", str(out_dir)]),
        ("Anderson-Darling normality audit", "analyze_ad_norm.py",
         ["--data", str(data_path), "--out", str(out_dir)]),
        ("Logistic discrimination modeling & audit", "discrimination.py",
         ["--data", str(data_path), "--out", str(out_dir)]),
        ("Treatment plan characteristics summary", "plan_characteristics.py",
         ["--data", str(data_path), "--out", str(out_dir)]),
        ("Summary figure generation", "figures.py",
         ["--data", str(data_path), "--results", str(out_dir), "--out", str(fig_dir)]),
    ]

    total_start = time.perf_counter()
    print("=" * 70)
    print("  Helical Tomotherapy PSQA: Reproducing Study Analysis & Figures")
    print(f"  Input data: {data_path.name}")
    print(f"  Output directory: {out_dir}")
    print("=" * 70)

    for i, (desc, script_name, options) in enumerate(pipeline, start=1):
        print(f"\n[{i}/{len(pipeline)}] {desc} ({script_name})...", flush=True)
        step_start = time.perf_counter()
        script_path = ROOT / "analysis" / script_name
        subprocess.run([sys.executable, str(script_path), *options], check=True)
        print(f"  Done in {time.perf_counter() - step_start:.1f}s")

    elapsed = time.perf_counter() - total_start
    print("\n" + "=" * 70)
    print(f"Pipeline successfully completed in {elapsed:.1f}s.")
    print(f"All generated tables and figures are located in:\n  {out_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
