

# Interactive
srun --mem=5G -c 15 -J admixpop_WP1 --account=admixpop -t 8:00:00 --pty bash

# Get status of slurm jobs with full jobnames
squeue -u $USER -o "%.10i %.20P %.40j %.8T %.10M %.10l %.10m %R"





