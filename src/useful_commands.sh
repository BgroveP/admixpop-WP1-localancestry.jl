

# Interactive
srun --mem=5G -c 15 -J admixpop_WP1 --account=admixpop -t 8:00:00 --pty bash

# Get status of slurm jobs with full jobnames
squeue -u $USER -o "%.10i %.20P %.40j %.8T %.10M %.10l %.10m %R"


# test

export CMAKE_PREFIX_PATH="$CONDA_PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export LIBRARY_PATH="$CONDA_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CPATH="$CONDA_PREFIX/include${CPATH:+:$CPATH}"

SparsePainter/SparsePainter \
    -reffile data/replicate_1/ancestral_array_10_10.vcf.gz \
    -targetfile data/replicate_1/admixed_array_10_10.vcf.gz \
    -mapfile data/sparsepainter_array.tsv \
    -popfile data/replicate_1/global_ancestries_sparsepainter.tsv \
    -namefile data/replicate_1/individuals_admixed_gen10_10.txt \
    -out results/localancestry/sparsepainter/1_array_10_10 \
    -diff_lambda \
    -prob \
    -probstore raw \
    -ncores 10

julia -t 10 --project=. src/julia/accuracy_sparsepainter.jl \
    results/localancestry/sparsepainter/1_array_10_10_prob.txt.gz \
    data/replicate_1/local_ancestries.csv \
    results/accuracy/sparsepainter/1_array_10_10.csv

