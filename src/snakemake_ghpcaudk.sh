#!/bin/bash

mkdir -p logs/snakemake
nohup pixi run snakemake --profile profiles/ghpcaudk --rerun-incomplete > logs/snakemake/$(date +%Y%m%d_%H%M%S).txt 2>&1 &


