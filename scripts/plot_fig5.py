#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def taylor_law(ca: float) -> float:
    c23 = ca ** (2.0 / 3.0)
    return (1.34 * c23) / (1.0 + 1.34 * 2.5 * c23)


def main():
    parser = argparse.ArgumentParser(description="Plot clean Taylor/Bretherton benchmark (Fig. 5 style).")
    parser.add_argument("--csv", default="intermediate/clean_sweep.csv")
    parser.add_argument("--out", default="intermediate/clean_benchmark_fig5.png")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f"CSV not found: {csv_path}")

    ca_vals = []
    h_vals = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ca_vals.append(float(row["Ca"]))
            h_vals.append(float(row["h_over_R"]))

    if not ca_vals:
        raise SystemExit("No data in CSV")

    ca_min = min(ca_vals)
    ca_max = max(ca_vals)

    # Smooth curves for Taylor and upper limit
    npts = 200
    ca_curve = [ca_min * (ca_max / ca_min) ** (i / (npts - 1)) for i in range(npts)]
    taylor_curve = [taylor_law(ca) for ca in ca_curve]
    upper_curve = [(4 ** (2.0 / 3.0)) * h for h in taylor_curve]

    fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=160)
    ax.loglog(ca_vals, h_vals, "o", label="Clean (Basilisk)")
    ax.loglog(ca_curve, taylor_curve, "-", label="Taylor (Aussillous & Quéré)")
    ax.loglog(ca_curve, upper_curve, "-.", label=r"$4^{2/3}$ Taylor")

    ax.set_xlabel("Ca")
    ax.set_ylabel(r"$h_\infty / R$")
    ax.grid(True, which="both", ls=":", lw=0.6)
    ax.legend()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
