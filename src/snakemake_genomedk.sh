#!/bin/bash
# Launch the full pipeline in the background on a SLURM cluster.
# Usage: bash src/snakemake_ghpcaudk.sh [profile]   (default: profiles/ghpcaudk)
# Only ONE Snakemake instance may run in this directory at a time.
set -euo pipefail

PROFILE="${1:-profiles/genomedk}"
mkdir -p logs/snakemake logs/slurm
LOG="logs/snakemake/$(date +%Y%m%d_%H%M%S).txt"
nohup pixi run snakemake --profile "$PROFILE" --rerun-incomplete > "$LOG" 2>&1 &
echo "Started Snakemake with $PROFILE (PID $!). Log: $LOG"
