#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[INFO] plotting benchmark figures ..."
python3 scripts/plot_2dpaper_results.py --intermediate intermediate
echo "[INFO] plot done."
