#!/bin/bash
#SBATCH --mem=5G
#SBATCH -c 15
#SBATCH -J admixpop_WP1
#SBATCH --account=admixpop
#SBATCH --output=../../../docs/slurm/WP1_article.log 
#SBATCH --error=../../../docs/slurm/WP1_article.err 
#SBATCH -t 3-1:00:00 #Job max time

cd .. 
pixi run snakemake --cores $(nproc) --rerun-incomplete