#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=11:00:00
#SBATCH --cpus-per-task=24
#SBATCH --output=tmp/rfmix_timetest.log

time rfmix/rfmix --n-threads=24 \
  -f data/replicate_1/admixed_sequence_6_10.vcf.gz \
  -r data/replicate_1/ancestral_sequence_6_10.vcf.gz \
  -m data/replicate_1/global_ancestries_rfmix.tsv \
  -g data/rfmix_sequence.tsv \
  -o tmp/rfmix_timetest \
  --chromosome=chr22