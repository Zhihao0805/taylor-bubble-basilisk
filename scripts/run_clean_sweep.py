#!/usr/bin/env python3
import argparse
import csv
import math
import os
import subprocess
import sys
from pathlib import Path

CA_VALUES = [0.0015, 0.004, 0.016, 0.044, 0.097]


def taylor_law(ca: float) -> float:
    c23 = ca ** (2.0 / 3.0)
    return (1.34 * c23) / (1.0 + 1.34 * 2.5 * c23)


def parse_result(intermediate: Path) -> dict:
    final_path = intermediate / "final.csv"
    summary_path = intermediate / "summary.csv"

    if final_path.exists():
        with open(final_path, newline="") as f:
            row = next(csv.DictReader(f), None)
        if not row:
            raise RuntimeError(f"final.csv is empty: {final_path}")
        return {
            "h_over_R": float(row["h_inf_over_R"]),
            "h_taylor_over_R": float(row["h_Taylor_over_R"]),
            "rel_error": float(row["rel_error"]),
            "Ub_mean": float(row["Ub_mean"]),
            "Vgas_mean": float(row["Vgas_mean"]),
        }

    if summary_path.exists():
        with open(summary_path, newline="") as f:
            rows = list(csv.DictReader(f))
        if not rows:
            raise RuntimeError(f"summary.csv is empty: {summary_path}")
        row = rows[-1]
        if "h_over_R" in row:
            return {
                "h_over_R": float(row["h_over_R"]),
                "h_taylor_over_R": float(row["h_taylor_over_R"]),
                "rel_error": float(row["rel_error"]),
            }
        if "h_inf_over_R" in row:
            h = float(row["h_inf_over_R"])
            ca = float(row["Ca"])
            h_t = taylor_law(ca)
            return {
                "h_over_R": h,
                "h_taylor_over_R": h_t,
                "rel_error": abs(h - h_t) / max(1e-30, h_t),
                "Ub_mean": float(row["Ub"]),
                "Vgas_mean": float(row["Vgas"]),
            }
        raise RuntimeError(f"summary.csv format not recognized: {summary_path}")

    raise RuntimeError(f"No final.csv or summary.csv found in {intermediate}")


def run_case(run_exe: Path, ca: float, args) -> dict:
    intermediate = Path(args.intermediate)
    intermediate.mkdir(parents=True, exist_ok=True)

    maxlevel = args.maxlevel_lowca if ca < args.lowca_threshold else args.maxlevel
    cmd = [
        str(run_exe),
        str(maxlevel),
        str(args.re),
        str(ca),
        str(args.lz_over_r),
        str(args.tmax),
    ]

    log_out = intermediate / f"run_Ca_{ca:.4g}.log"
    log_err = intermediate / f"run_Ca_{ca:.4g}.err"

    with open(log_out, "w") as f_out, open(log_err, "w") as f_err:
        subprocess.run(cmd, stdout=f_out, stderr=f_err, check=True)

    parsed = parse_result(intermediate)
    h_over_r = parsed["h_over_R"]
    h_taylor = parsed["h_taylor_over_R"]
    rel_error = parsed["rel_error"]

    return {
        "Ca": ca,
        "MAXLEVEL": maxlevel,
        "h_over_R": h_over_r,
        "h_taylor_over_R": h_taylor,
        "rel_error": rel_error,
        "Ub_mean": parsed.get("Ub_mean", float("nan")),
        "Vgas_mean": parsed.get("Vgas_mean", float("nan")),
    }


def main():
    parser = argparse.ArgumentParser(description="Run clean Taylor/Bretherton Ca sweep.")
    parser.add_argument("--run", default="./run", help="Path to Basilisk executable (default: ./run)")
    parser.add_argument("--maxlevel", type=int, default=12)
    parser.add_argument("--maxlevel-lowca", type=int, default=13)
    parser.add_argument("--lowca-threshold", type=float, default=0.004)
    parser.add_argument("--ca-values", default="0.0015,0.004,0.016,0.044,0.097")
    parser.add_argument("--re", type=float, default=1.0)
    parser.add_argument("--lz-over-r", type=float, default=40.0)
    parser.add_argument("--tmax", type=float, default=2.0)
    parser.add_argument("--intermediate", default="intermediate")
    parser.add_argument("--csv", default="intermediate/clean_sweep.csv")
    args = parser.parse_args()

    run_exe = Path(args.run)
    if not run_exe.is_absolute():
        run_exe = (Path.cwd() / run_exe).resolve()
    if not run_exe.exists():
        print(f"ERROR: run executable not found at {run_exe}", file=sys.stderr)
        print("Build it first, e.g.: qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run", file=sys.stderr)
        return 1

    ca_values = [float(x.strip()) for x in args.ca_values.split(",") if x.strip()]
    results = []
    for ca in ca_values:
        result = run_case(run_exe, ca, args)
        # Override Taylor/rel_error if needed to ensure consistent formula
        h_taylor = taylor_law(ca)
        rel_error = abs(result["h_over_R"] - h_taylor) / max(1e-30, h_taylor)
        result["h_taylor_over_R"] = h_taylor
        result["rel_error"] = rel_error
        results.append(result)
        print(
            f"Ca={ca:.4g} MAXLEVEL={result['MAXLEVEL']} "
            f"h/R={result['h_over_R']:.6g} Taylor={h_taylor:.6g} "
            f"rel_error={rel_error:.4g} Ub_mean={result['Ub_mean']:.6g}"
        )

    out_csv = Path(args.csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Ca",
                "MAXLEVEL",
                "h_over_R",
                "h_taylor_over_R",
                "rel_error",
                "Ub_mean",
                "Vgas_mean",
            ],
        )
        writer.writeheader()
        writer.writerows(results)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
