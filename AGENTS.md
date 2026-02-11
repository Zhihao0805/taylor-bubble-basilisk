# Project rules (Basilisk Taylor bubble)

## Build & run
- Compiler: qcc (Basilisk)
- Primary case file: taylor_clean_benchmark_axi_steady.c (or your main .c)
- Run script: ./run MAXLEVEL Re Ca Lz_over_R TMAX
- Always keep output columns unchanged unless explicitly requested.

## Constraints
- Must support embedded boundary (embed.h) + tree AMR without crashing.
- If using AMR: after adapt, rebuild cs/fs from global vertex scalar phi.
- Prefer portable Basilisk includes (adapt.h vs adapt_wavelet.h).
- Do not add external dependencies beyond standard Linux utils.

## Deliverables
- A clean-case benchmark matching docs/clean_benchmark.txt:
  geometry, nondimensionalization, BCs, diagnostics, convergence checks.
- Reproducible run instructions in README snippet.
