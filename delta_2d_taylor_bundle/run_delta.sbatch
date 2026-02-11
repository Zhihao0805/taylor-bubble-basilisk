#!/bin/bash
#SBATCH --job-name=taylor2d-paper
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

# If your Delta environment needs module setup, add it here.
# module purge
# module load gcc python

bash delta_2d_taylor_bundle/01_preprocess_build.sh
bash delta_2d_taylor_bundle/02_run_campaigns.sh
bash delta_2d_taylor_bundle/03_postprocess_plot.sh
bash delta_2d_taylor_bundle/04_check_outputs.sh

echo "[INFO] completed on Delta. If results are good, run:"
echo "bash delta_2d_taylor_bundle/05_autogit.sh \"Delta run: 2D Taylor bubble campaigns\""
