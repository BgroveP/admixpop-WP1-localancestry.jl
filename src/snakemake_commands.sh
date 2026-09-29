

# Interactive
srun --mem=5G -c 15 -J admixpop_WP1 --account=admixpop -t 8:00:00 --pty bash

# Run certain rule
pixi run snakemake --profile profiles/slurm --allowed-rules estimate_gnomix all estimate_localancestryjl
