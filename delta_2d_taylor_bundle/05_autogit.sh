#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

MSG="${1:-Delta bundle: 2D Taylor bubble preprocess+postprocess scripts}"

echo "[INFO] running autopush with message: $MSG"
./scripts/autopush.sh "$MSG"
