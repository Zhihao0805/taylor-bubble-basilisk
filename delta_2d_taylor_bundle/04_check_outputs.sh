#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

required=(
  intermediate/grid_convergence.csv
  intermediate/tol_sensitivity.csv
  intermediate/ca_sweep.csv
  intermediate/grid_convergence.png
  intermediate/tol_sensitivity.png
  intermediate/ca_sweep_vs_theory.png
)

echo "[INFO] checking required outputs ..."
for f in "${required[@]}"; do
  if [ ! -f "$f" ]; then
    echo "[FAIL] missing: $f"
    exit 1
  fi
  echo "[OK] $f"
done

mp4="$(ls intermediate/movie_*_centered.mp4 2>/dev/null | head -n 1 || true)"
if [ -z "$mp4" ]; then
  mp4="$(ls intermediate/*_centered.mp4 2>/dev/null | head -n 1 || true)"
fi
if [ -z "$mp4" ]; then
  echo "[FAIL] centered mp4 not found under intermediate/"
  exit 1
fi

if command -v ffprobe >/dev/null 2>&1; then
  echo "[INFO] ffprobe: $mp4"
  ffprobe -v error -count_frames -select_streams v:0 \
    -show_entries stream=nb_read_frames,duration,width,height \
    -of default=nokey=1:noprint_wrappers=1 "$mp4"
else
  echo "[WARN] ffprobe not found; skip mp4 metadata check."
fi

echo "[INFO] git status (should not stage runtime artifacts):"
git status --short
echo "[PASS] output check done."
