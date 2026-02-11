#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def load_rows(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def fval(r, key):
    return float(r[key])


def plot_grid(grid_csv: Path, out_png: Path):
    rows = load_rows(grid_csv)
    rows.sort(key=lambda r: float(r["H_over_Delta"]))
    x = [fval(r, "H_over_Delta") for r in rows]
    y = [fval(r, "Ud_over_Uf_mean_last30") for r in rows]
    y_a = [fval(r, "Ud_over_Uf_aussillous") for r in rows]

    fig, ax = plt.subplots(figsize=(6.2, 4.2), dpi=160)
    ax.plot(x, y, "o-", label="Simulation")
    if y_a:
        ax.plot(x, [y_a[-1]] * len(x), "--", label="Aussillous (Eq. 8-9)")
    ax.set_xlabel(r"$H/\Delta$")
    ax.set_ylabel(r"$U_d/U_f$")
    ax.set_title("Grid Convergence (Ca=0.001, Re=0.1)")
    ax.grid(True, ls=":", lw=0.6)
    ax.legend()
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)


def plot_tol(ts_csv: Path, out_png: Path):
    rows = load_rows(ts_csv)
    groups = defaultdict(list)
    for r in rows:
        groups[r["tol"]].append(r)

    fig, ax = plt.subplots(figsize=(6.2, 4.2), dpi=160)
    for tol, rr in sorted(groups.items(), key=lambda kv: float(kv[0])):
        rr.sort(key=lambda r: float(r["t"]))
        t = [fval(r, "t") for r in rr]
        y = [fval(r, "Ud_over_Uf") for r in rr]
        ax.plot(t, y, label=f"TOL={tol}")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$U_d/U_f$")
    ax.set_title("Poisson Tolerance Sensitivity")
    ax.grid(True, ls=":", lw=0.6)
    ax.legend()
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)


def plot_ca(ca_csv: Path, out_png: Path):
    rows = load_rows(ca_csv)
    rows.sort(key=lambda r: float(r["Ca"]))
    ca = [fval(r, "Ca") for r in rows]
    y_sim = [fval(r, "Ud_over_Uf_sim") for r in rows]
    y_b = [fval(r, "Ud_over_Uf_bretherton") for r in rows]
    y_a = [fval(r, "Ud_over_Uf_aussillous") for r in rows]

    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=160)
    ax.semilogx(ca, y_sim, "o", label="Simulation")
    ax.semilogx(ca, y_b, "-", label="Bretherton (Eq. 7-9)")
    ax.semilogx(ca, y_a, "--", label="Aussillous (Eq. 8-9)")
    ax.set_xlabel("Ca")
    ax.set_ylabel(r"$U_d/U_f$")
    ax.set_title("Ca Sweep vs Theory (Re=0.1)")
    ax.grid(True, which="both", ls=":", lw=0.6)
    ax.legend()
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)


def main():
    ap = argparse.ArgumentParser(description="Plot Section-4 2D Taylor bubble benchmark results.")
    ap.add_argument("--intermediate", default="intermediate")
    ap.add_argument("--grid-csv", default=None)
    ap.add_argument("--tol-ts-csv", default=None)
    ap.add_argument("--ca-csv", default=None)
    ap.add_argument("--grid-out", default=None)
    ap.add_argument("--tol-out", default=None)
    ap.add_argument("--ca-out", default=None)
    args = ap.parse_args()

    base = Path(args.intermediate)
    grid_csv = Path(args.grid_csv) if args.grid_csv else base / "grid_convergence.csv"
    tol_csv = Path(args.tol_ts_csv) if args.tol_ts_csv else base / "tol_timeseries.csv"
    ca_csv = Path(args.ca_csv) if args.ca_csv else base / "ca_sweep.csv"
    grid_out = Path(args.grid_out) if args.grid_out else base / "grid_convergence.png"
    tol_out = Path(args.tol_out) if args.tol_out else base / "tol_sensitivity.png"
    ca_out = Path(args.ca_out) if args.ca_out else base / "ca_sweep_vs_theory.png"

    plot_grid(grid_csv, grid_out)
    plot_tol(tol_csv, tol_out)
    plot_ca(ca_csv, ca_out)

    print(f"Saved: {grid_out}")
    print(f"Saved: {tol_out}")
    print(f"Saved: {ca_out}")


if __name__ == "__main__":
    main()
