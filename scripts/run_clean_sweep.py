#!/usr/bin/env python3
import argparse
import csv
import math
import os
import subprocess
import sys
from pathlib import Path

CA_VALUES = [0.0015, 0.004, 0.008, 0.016, 0.044, 0.097]


def taylor_law(ca: float) -> float:
    c23 = ca ** (2.0 / 3.0)
    return (1.34 * c23) / (1.0 + 1.34 * 2.5 * c23)


def run_case(run_exe: Path, ca: float, args) -> dict:
    intermediate = Path(args.intermediate)
    intermediate.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(run_exe),
        str(args.maxlevel),
        str(args.re),
        str(ca),
        str(args.lz_over_r),
        str(args.tmax),
    ]

    log_out = intermediate / f"run_Ca_{ca:.4g}.log"
    log_err = intermediate / f"run_Ca_{ca:.4g}.err"

    with open(log_out, "w") as f_out, open(log_err, "w") as f_err:
        subprocess.run(cmd, stdout=f_out, stderr=f_err, check=True)

    summary_path = intermediate / "summary.csv"
    if not summary_path.exists():
        raise RuntimeError(f"summary.csv not found after run: {summary_path}")

    with open(summary_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        if not rows:
            raise RuntimeError("summary.csv is empty")
        row = rows[0]

    h_over_r = float(row["h_over_R"])
    h_taylor = float(row["h_taylor_over_R"])
    rel_error = float(row["rel_error"])

    return {
        "Ca": ca,
        "h_over_R": h_over_r,
        "h_taylor_over_R": h_taylor,
        "rel_error": rel_error,
    }


def main():
    parser = argparse.ArgumentParser(description="Run clean Taylor/Bretherton Ca sweep.")
    parser.add_argument("--run", default="./run", help="Path to Basilisk executable (default: ./run)")
    parser.add_argument("--maxlevel", type=int, default=9)
    parser.add_argument("--re", type=float, default=1.0)
    parser.add_argument("--lz-over-r", type=float, default=40.0)
    parser.add_argument("--tmax", type=float, default=2.0)
    parser.add_argument("--intermediate", default="intermediate")
    parser.add_argument("--csv", default="intermediate/clean_sweep.csv")
    args = parser.parse_args()

    run_exe = Path(args.run)
    if not run_exe.exists():
        print(f"ERROR: run executable not found at {run_exe}", file=sys.stderr)
        print("Build it first, e.g.: qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run", file=sys.stderr)
        return 1

    results = []
    for ca in CA_VALUES:
        result = run_case(run_exe, ca, args)
        # Override Taylor/rel_error if needed to ensure consistent formula
        h_taylor = taylor_law(ca)
        rel_error = abs(result["h_over_R"] - h_taylor) / max(1e-30, h_taylor)
        result["h_taylor_over_R"] = h_taylor
        result["rel_error"] = rel_error
        results.append(result)
        print(f"Ca={ca:.4g} h/R={result['h_over_R']:.6g} Taylor={h_taylor:.6g} rel_error={rel_error:.4g}")

    out_csv = Path(args.csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Ca", "h_over_R", "h_taylor_over_R", "rel_error"])
        writer.writeheader()
        writer.writerows(results)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
