#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from pathlib import Path


def eq7_two_t_over_h(ca: float) -> float:
    return 0.643 * (3.0 * ca) ** (2.0 / 3.0)


def eq8_two_t_over_h(ca: float) -> float:
    num = eq7_two_t_over_h(ca)
    return num / (1.0 + 0.643 * 2.5 * (3.0 * ca) ** (2.0 / 3.0))


def eq9_uduf(two_t_over_h: float) -> float:
    denom = max(1e-12, 1.0 - two_t_over_h)
    return 1.0 / denom


def slug(v: str) -> str:
    return v.replace(".", "p").replace("-", "m").replace("+", "")


def run_case(run_exe: Path, intermediate: Path, lev: int, re: float, ca: float, tmax: float,
             tol: float, tag: str, make_movie: int, eps_u: float, eps_f: float, write_dumps: int):
    cmd = [
        str(run_exe),
        str(lev),
        str(re),
        str(ca),
        str(tmax),
        f"{tol:.1e}",
        tag,
        str(make_movie),
        str(eps_u),
        str(eps_f),
        str(write_dumps),
    ]

    out_log = intermediate / f"{tag}.out"
    err_log = intermediate / f"{tag}.err"
    with out_log.open("w") as fout, err_log.open("w") as ferr:
        subprocess.run(cmd, stdout=fout, stderr=ferr, check=True)

    final_csv = intermediate / f"{tag}_final.csv"
    if not final_csv.exists():
        raise RuntimeError(f"missing final csv: {final_csv}")
    with final_csv.open(newline="") as f:
        row = next(csv.DictReader(f), None)
    if row is None:
        raise RuntimeError(f"empty final csv: {final_csv}")

    return row


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Section-4 2D Taylor bubble benchmark campaigns.")
    parser.add_argument("--run", default="./run2d", help="solver executable path")
    parser.add_argument("--build-cmd", default="qcc -O2 -Wall -disable-dimensions taylor_benchmark_2Dpaper.c -lm -o run2d")
    parser.add_argument("--intermediate", default="intermediate")
    parser.add_argument("--re", type=float, default=0.1)
    parser.add_argument("--grid-ca", type=float, default=0.001)
    parser.add_argument("--grid-levs", default="8,9,10,11,12")
    parser.add_argument("--grid-tmax", type=float, default=0.5)
    parser.add_argument("--tol-values", default="1e-3,1e-5,1e-7")
    parser.add_argument("--tol-lev", type=int, default=11)
    parser.add_argument("--tol-ca", type=float, default=0.001)
    parser.add_argument("--tol-tmax", type=float, default=0.8)
    parser.add_argument("--ca-values", default="0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1")
    parser.add_argument("--ca-lev", type=int, default=11)
    parser.add_argument("--ca-lev-low", type=int, default=12)
    parser.add_argument("--ca-low-threshold", type=float, default=0.0006)
    parser.add_argument("--ca-tmax", type=float, default=0.5)
    parser.add_argument("--eps-u", type=float, default=7e-3)
    parser.add_argument("--eps-f", type=float, default=1e-1)
    parser.add_argument("--tol-default", type=float, default=1e-5)
    parser.add_argument("--movie-ca", type=float, default=0.01)
    parser.add_argument("--movie-lev", type=int, default=11)
    parser.add_argument("--movie-tmax", type=float, default=0.8)
    args = parser.parse_args()

    root = Path.cwd()
    intermediate = Path(args.intermediate)
    intermediate.mkdir(parents=True, exist_ok=True)

    run_exe = Path(args.run)
    if not run_exe.is_absolute():
        run_exe = (root / run_exe).resolve()
    if not run_exe.exists():
        subprocess.run(args.build_cmd, shell=True, check=True)
        if not run_exe.exists():
            raise RuntimeError(f"build finished but executable not found: {run_exe}")

    re = args.re

    # Campaign 1: grid convergence
    levs = [int(v.strip()) for v in args.grid_levs.split(",") if v.strip()]
    grid_rows = []
    for lev in levs:
        tag = f"grid_lev{lev}_re{slug(f'{re:g}')}_ca{slug(f'{args.grid_ca:g}')}_tol{slug(f'{args.tol_default:.1e}')}"
        row = run_case(run_exe, intermediate, lev, re, args.grid_ca, args.grid_tmax, args.tol_default,
                       tag, make_movie=0, eps_u=args.eps_u, eps_f=args.eps_f, write_dumps=0)
        row["campaign"] = "grid_convergence"
        row["H_over_Delta"] = str((2 ** lev) / 10.0)  # H=1, L0=10
        grid_rows.append(row)
        print(f"[grid] LEV={lev} H/Delta={(2 ** lev)/10.0:.2f} Ud/Uf={float(row['Ud_over_Uf_mean_last30']):.6g}")

    grid_csv = intermediate / "grid_convergence.csv"
    grid_fields = [
        "campaign", "tag", "Re", "Ca", "LEV", "TOL", "H_over_Delta",
        "Ud_over_Uf_mean_last30", "Ud_over_Uf_aussillous", "Ud_over_Uf_bretherton",
        "rel_err_aussillous", "rel_err_bretherton", "n_steady_samples"
    ]
    with grid_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=grid_fields)
        w.writeheader()
        for r in grid_rows:
            w.writerow({k: r.get(k, "") for k in grid_fields})

    # Campaign 2: Poisson tolerance sensitivity
    tol_values = [float(v.strip()) for v in args.tol_values.split(",") if v.strip()]
    tol_rows = []
    tol_ts_csv = intermediate / "tol_timeseries.csv"
    ts_fields = None
    with tol_ts_csv.open("w", newline="") as fts:
        ts_writer = None
        for tol in tol_values:
            tag = f"tol_lev{args.tol_lev}_re{slug(f'{re:g}')}_ca{slug(f'{args.tol_ca:g}')}_tol{slug(f'{tol:.1e}')}"
            row = run_case(run_exe, intermediate, args.tol_lev, re, args.tol_ca, args.tol_tmax, tol,
                           tag, make_movie=0, eps_u=args.eps_u, eps_f=args.eps_f, write_dumps=0)
            row["campaign"] = "tol_sensitivity"
            tol_rows.append(row)
            print(f"[tol] TOL={tol:.1e} Ud/Uf={float(row['Ud_over_Uf_mean_last30']):.6g}")

            # merge timeseries for plotting Ud/Uf(t)
            ts_path = intermediate / f"{tag}_timeseries.csv"
            with ts_path.open(newline="") as fcase:
                reader = csv.DictReader(fcase)
                if ts_writer is None:
                    ts_fields = ["tag", "tol"] + reader.fieldnames
                    ts_writer = csv.DictWriter(fts, fieldnames=ts_fields)
                    ts_writer.writeheader()
                for rr in reader:
                    rr_out = {"tag": tag, "tol": f"{tol:.1e}"}
                    rr_out.update(rr)
                    ts_writer.writerow(rr_out)

    tol_csv = intermediate / "tol_sensitivity.csv"
    tol_fields = [
        "campaign", "tag", "Re", "Ca", "LEV", "TOL",
        "Ud_over_Uf_mean_last30", "Ud_over_Uf_aussillous", "Ud_over_Uf_bretherton",
        "rel_err_aussillous", "rel_err_bretherton", "n_steady_samples"
    ]
    with tol_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tol_fields)
        w.writeheader()
        for r in tol_rows:
            w.writerow({k: r.get(k, "") for k in tol_fields})

    # Campaign 3: Ca sweep
    ca_values = [float(v.strip()) for v in args.ca_values.split(",") if v.strip()]
    ca_rows = []
    for ca in ca_values:
        lev = args.ca_lev_low if ca <= args.ca_low_threshold else args.ca_lev
        tag = f"ca_lev{lev}_re{slug(f'{re:g}')}_ca{slug(f'{ca:g}')}_tol{slug(f'{args.tol_default:.1e}')}"
        row = run_case(run_exe, intermediate, lev, re, ca, args.ca_tmax, args.tol_default,
                       tag, make_movie=0, eps_u=args.eps_u, eps_f=args.eps_f, write_dumps=0)
        two_t_b = eq7_two_t_over_h(ca)
        two_t_a = eq8_two_t_over_h(ca)
        uduf_b = eq9_uduf(two_t_b)
        uduf_a = eq9_uduf(two_t_a)
        uduf_sim = float(row["Ud_over_Uf_mean_last30"])
        rel_a = abs(uduf_sim - uduf_a) / max(1e-12, uduf_a)
        rel_b = abs(uduf_sim - uduf_b) / max(1e-12, uduf_b)
        out = {
            "campaign": "ca_sweep",
            "tag": tag,
            "Re": re,
            "Ca": ca,
            "LEV": lev,
            "TOL": args.tol_default,
            "Ud_over_Uf_sim": uduf_sim,
            "Ud_over_Uf_bretherton": uduf_b,
            "Ud_over_Uf_aussillous": uduf_a,
            "rel_err_bretherton": rel_b,
            "rel_err_aussillous": rel_a,
            "n_steady_samples": row.get("n_steady_samples", ""),
        }
        ca_rows.append(out)
        print(f"[ca] Ca={ca:.4g} LEV={lev} Ud/Uf={uduf_sim:.6g} err_a={rel_a:.4g}")

    ca_csv = intermediate / "ca_sweep.csv"
    ca_fields = [
        "campaign", "tag", "Re", "Ca", "LEV", "TOL",
        "Ud_over_Uf_sim", "Ud_over_Uf_bretherton", "Ud_over_Uf_aussillous",
        "rel_err_bretherton", "rel_err_aussillous", "n_steady_samples"
    ]
    with ca_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ca_fields)
        w.writeheader()
        w.writerows(ca_rows)

    # Representative centered movie case
    movie_tag = f"movie_lev{args.movie_lev}_re{slug(f'{re:g}')}_ca{slug(f'{args.movie_ca:g}')}_tol{slug(f'{args.tol_default:.1e}')}"
    run_case(run_exe, intermediate, args.movie_lev, re, args.movie_ca, args.movie_tmax, args.tol_default,
             movie_tag, make_movie=1, eps_u=args.eps_u, eps_f=args.eps_f, write_dumps=0)
    print(f"[movie] generated: {intermediate / (movie_tag + '_centered.mp4')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
