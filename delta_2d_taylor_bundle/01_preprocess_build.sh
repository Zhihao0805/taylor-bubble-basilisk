#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if ! command -v qcc >/dev/null 2>&1; then
  echo "[ERROR] qcc not found in PATH. Please load Basilisk environment first."
  exit 1
fi

TS="$(date +%Y%m%d-%H%M%S)"
BKP="backups/${TS}"
mkdir -p "$BKP"

if [ -d intermediate ]; then
  mv intermediate "$BKP/intermediate"
fi
mkdir -p intermediate

echo "[INFO] backup_dir=$BKP"
echo "[INFO] building run2d ..."
qcc -O2 -Wall -disable-dimensions delta_2d_taylor_bundle/taylor_benchmark_2Dpaper_bundle.c -lm -o run2d
echo "[INFO] build done: $ROOT_DIR/run2d"
