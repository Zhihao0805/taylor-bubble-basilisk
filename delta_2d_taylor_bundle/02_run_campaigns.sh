#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

source delta_2d_taylor_bundle/params_delta.env

if [ ! -x ./run2d ]; then
  echo "[ERROR] run2d not found. Run: bash delta_2d_taylor_bundle/01_preprocess_build.sh"
  exit 1
fi

echo "[INFO] running campaigns ..."
PYTHONUNBUFFERED=1 python3 scripts/run_2dpaper_campaigns.py \
  --run ./run2d \
  --re "$RE" \
  --grid-ca "$GRID_CA" \
  --grid-levs "$GRID_LEVS" \
  --grid-tmax "$GRID_TMAX" \
  --tol-values "$TOL_VALUES" \
  --tol-lev "$TOL_LEV" \
  --tol-ca "$TOL_CA" \
  --tol-tmax "$TOL_TMAX" \
  --ca-values "$CA_VALUES" \
  --ca-lev "$CA_LEV" \
  --ca-lev-low "$CA_LEV_LOW" \
  --ca-low-threshold "$CA_LOW_THRESHOLD" \
  --ca-tmax "$CA_TMAX" \
  --tol-default "$TOL_DEFAULT" \
  --eps-u "$EPS_U" \
  --eps-f "$EPS_F" \
  --movie-ca "$MOVIE_CA" \
  --movie-lev "$MOVIE_LEV" \
  --movie-tmax "$MOVIE_TMAX"

echo "[INFO] campaigns done."
