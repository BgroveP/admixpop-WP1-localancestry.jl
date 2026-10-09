#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=3-11:00:00
#SBATCH --cpus-per-task=24
#SBATCH --output=tmp/sparsepainter_timetest.log
#SBATCH -p zen4

eval "$(pixi shell-hook)"

# Expose LAPACK and BLAS
export CMAKE_PREFIX_PATH="$CONDA_PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export LIBRARY_PATH="$CONDA_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CPATH="$CONDA_PREFIX/include${CPATH:+:$CPATH}"

time SparsePainter/SparsePainter \
    -reffile data/replicate_1/ancestral_sequence_6_10.vcf.gz \
    -targetfile data/replicate_1/admixed_sequence_6_10.vcf.gz \
    -mapfile data/sparsepainter_sequence.tsv \
    -popfile data/replicate_1/global_ancestries_sparsepainter_6.tsv \
    -namefile data/replicate_1/individuals_admixed_gen6_10.txt \
    -out results/localancestry/sparsepainter/1_sequence_6_6 \
    -prob \
    -probstore raw \
    -ncores 24