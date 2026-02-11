# Clean Benchmark Log (2026-02-10)

## Reference benchmark source
- Source: `Clean_benchmark.pdf` (paper: *Effects of surfactant on liquid film thickness in the Bretherton problem*).
- Clean-case benchmark statements used:
  - Axisymmetric tube geometry; axial length requires long domain, with `Lz/R` up to `40` for steady motion.
  - Clean-case target trend against Taylor law:
    - `h/R = 1.34*Ca^(2/3) / (1 + 2.5*1.34*Ca^(2/3))`
  - Base nondimensional setup references `Re = 1` and `Lz/R = 40`.
  - Grid guidance in the paper: radial resolution >= 64, and around 128 for very low `Ca` (`Ca < 0.004`).

## Repo markdown instructions discovered
Files checked:
- `AGENTS.md`
- `terminal.md`

Exact workflow commands extracted from markdown:
```bash
qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run
./run MAXLEVEL Re Ca Lz_over_R TMAX
./run 9 1 0.016 40 2 | tee intermediate/run.log
python3 scripts/run_clean_sweep.py
python3 scripts/plot_fig5.py
```

## Backup safety actions (no deletions)
Timestamped backups created and prior outputs moved before new runs:
- `backups/20260210-214052/`
- `backups/20260210-214338/`
- `backups/20260210-215019/`
- `backups/20260210-215520/`
- `backups/20260210-215713/`

What was transferred:
- prior `intermediate/` contents (dumps, logs, csv, old artifacts)
- intermediate outputs from aborted/partial runs

Clean-start policy used each time:
- move existing `intermediate/` into `backups/<timestamp>/intermediate/`
- recreate empty `intermediate/` before next run

## Parameter set used for CLEAN benchmark workflow
- Geometry/domain: axisymmetric tube, `Lz/R = 40`
- Nondimensional groups: `Re = 1.0`
- Ca list (5-point clean sweep):
  - `0.0015, 0.004, 0.016, 0.044, 0.097`
- Practical runtime sweep resolution used (best-effort low-cost verification):
  - `MAXLEVEL=7` for `Ca < 0.004`
  - `MAXLEVEL=6` otherwise
  - `TMAX=0.005` for sweep
- Full-duration single-case run for completion check:
  - `MAXLEVEL=6, Re=1, Ca=0.016, Lz/R=40, TMAX=2`

Steady criterion in current code:
- existing implementation keeps averaging diagnostics over the last 30% of simulation time (`t >= 0.7*t_end`) for summary output.

## Commands executed (compute -> postprocess -> visualize)
Build:
```bash
qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run
```

Short smoke test:
```bash
./run 9 1 0.016 40 0.2 | tee intermediate/testrun.log
```

5-point clean sweep:
```bash
PYTHONUNBUFFERED=1 python3 scripts/run_clean_sweep.py \
  --run ./run \
  --maxlevel 6 \
  --maxlevel-lowca 7 \
  --lowca-threshold 0.004 \
  --re 1.0 \
  --lz-over-r 40 \
  --tmax 0.005 \
  --intermediate intermediate \
  --csv intermediate/clean_sweep.csv \
  --ca-values 0.0015,0.004,0.016,0.044,0.097
```

Full-duration single case:
```bash
./run 6 1 0.016 40 2 | tee intermediate/full_run_Ca0016.log
```

Postprocess figure:
```bash
python3 scripts/plot_fig5.py --csv intermediate/clean_sweep.csv --out intermediate/clean_benchmark_fig5.png
```

MP4 verification:
```bash
ffprobe -v error -count_frames -select_streams v:0 \
  -show_entries stream=nb_read_frames,duration,width,height \
  -of default=nokey=1:noprint_wrappers=1 intermediate/clean_case.mp4
```

## Results
### 5-point sweep (`intermediate/clean_sweep.csv`)
| Ca | MAXLEVEL | h_sim/R | h_Taylor/R | rel_error |
|---:|---:|---:|---:|---:|
| 0.0015 | 7 | 0.218750 | 0.016821 | 12.0049 |
| 0.0040 | 6 | 0.062500 | 0.031137 | 1.0072 |
| 0.0160 | 6 | 0.062500 | 0.070161 | 0.1092 |
| 0.0440 | 6 | 0.062500 | 0.117818 | 0.4695 |
| 0.0970 | 6 | 0.062500 | 0.165702 | 0.6228 |

### Full-duration single case (`intermediate/final.csv`)
- `Ca=0.016`
- `h_inf/R=0.6875`
- `h_Taylor/R=0.0701606931899`
- `rel_error=8.79893397204`

## Visualization outputs
- MP4: `intermediate/clean_case.mp4`
- Figure: `intermediate/clean_benchmark_fig5.png`

`ffprobe` check for MP4:
- width: `600`
- height: `600`
- duration: `4.0 s`
- frame count: `100`

## Issues encountered and minimal fixes
1. Build include issue:
   - Error: `adapt.h` / `adapt_wavelet.h` not found in current Basilisk install.
   - Fix: removed non-portable adapt include from this case file and kept stable tree+embed path.

2. AMR stability issue:
   - Error: `embed-tree.h: Assertion 'coarse(cs)' failed` when AMR adaptation was enabled.
   - Fix: removed adaptive event in this run path (kept stable non-AMR execution).

3. Sweep script executable resolution bug:
   - Error: `FileNotFoundError: 'run'` when using `--run ./run`.
   - Fix: resolve run executable to absolute path before `subprocess.run()`.

4. Resolution initialization bug:
   - Issue: grid initialized with `N = 1 << MINLEVEL` causing coarse effective resolution.
   - Fix: switched to `N = 1 << MAXLEVEL` so CLI `MAXLEVEL` is honored.

5. Output/visualization workflow:
   - Added MP4 output event via `output_ppm(... file="intermediate/clean_case.mp4")`.
   - Added robust ignore entries to avoid staging runtime/binary artifacts.

## Notes
- This log captures a reproducible clean-case workflow with backup safety and artifact generation.
- Accuracy against paper-level clean benchmark remains resolution-sensitive; sweep above is a low-cost verification configuration.
