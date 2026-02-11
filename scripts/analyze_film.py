#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def taylor_law(ca: float) -> float:
    c23 = ca ** (2.0 / 3.0)
    return (1.34 * c23) / (1.0 + 1.34 * 2.5 * c23)


def load_summary(path: Path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})
    if not rows:
        raise RuntimeError(f"No data rows in {path}")
    return rows


def steady_stats(rows, frac=0.3):
    tmax = max(r["t"] for r in rows)
    tmin = tmax * (1.0 - frac)
    sel = [r for r in rows if r["t"] >= tmin]
    if not sel:
        sel = rows
    n = len(sel)
    h = sum(r["h_inf_over_R"] for r in sel) / n
    ub = sum(r["Ub"] for r in sel) / n
    vgas = sum(r["Vgas"] for r in sel) / n
    ca = sel[-1]["Ca"]
    return ca, h, ub, vgas, tmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="intermediate/summary.csv")
    ap.add_argument("--summary2", default=None, help="optional second-level summary.csv")
    args = ap.parse_args()

    rows = load_summary(Path(args.summary))
    ca, h, ub, vgas, tmax = steady_stats(rows)
    h_taylor = taylor_law(ca)
    rel_err = abs(h - h_taylor) / max(1e-30, h_taylor)

    print(f"steady window: last 30% of t (tmax={tmax:g})")
    print(f"Ca={ca:g}")
    print(f"h_inf/R={h:.6g}")
    print(f"Taylor={h_taylor:.6g}")
    print(f"rel_error={rel_err:.6g}")
    print(f"Ub_mean={ub:.6g}")
    print(f"Vgas_mean={vgas:.6g}")

    if args.summary2:
        rows2 = load_summary(Path(args.summary2))
        ca2, h2, _, _, tmax2 = steady_stats(rows2)
        if abs(ca2 - ca) > 1e-12:
            print(f"WARNING: Ca mismatch between summaries: {ca} vs {ca2}")
        rel_diff = abs(h2 - h) / max(1e-30, h2)
        print(f"level2 tmax={tmax2:g} h_inf/R={h2:.6g}")
        print(f"relative difference (level2 vs level1)={rel_diff:.6g}")


if __name__ == "__main__":
    main()
