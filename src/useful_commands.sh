

# Interactive
srun --mem=5G -c 15 -J admixpop_WP1 --account=admixpop -t 8:00:00 --pty bash

# Get status of slurm jobs with full jobnames
squeue -u $USER -o "%.10i %.20P %.40j %.8T %.10M %.10l %.10m %R"


#
# --- pass 1: initial panel estimation ---
java -jar flare/flare.jar \
    min-mac=0 min-maf=0.0 \
    ref=data/replicate_1/ancestral_array_6_10.vcf.gz \
    ref-panel=data/replicate_1/global_ancestries_flareclust_6.ssv \
    gt=data/replicate_1/admixed_array_6_10.vcf.gz \
    map=data/flare_array.ssv \
    panel-probs=true \
    out=results/localancestry/flareclust/1_array_6_6 \
    nthreads=10

# --- model file creation ---
python3 flare/create_model_file.py \
    3 \
    results/localancestry/flareclust/1_array_6_6.panels \
    results/localancestry/flareclust/1_array_6_6

# --- pass 2: ancestry estimation with fitted model ---
java -jar flare/flare.jar \
    min-mac=0 min-maf=0.0 \
    ref=data/replicate_1/ancestral_array_6_10.vcf.gz \
    ref-panel=data/replicate_1/global_ancestries_flareclust_6.ssv \
    gt=data/replicate_1/admixed_array_6_10.vcf.gz \
    map=data/flare_array.ssv \
    model=results/localancestry/flareclust/1_array_6_6.model \
    update-p=true \
    out=results/localancestry/flareclust/1_array_6_6 \
    nthreads=10
    
