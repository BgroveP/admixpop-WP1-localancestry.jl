
# LocalAncestry.jl: accuracy for a grid of maximum block sizes x IA thresholds.
# Uses the smallest panels (reference and target from generation 10 only).
# The script writes one file per (block size, threshold); the declared output is
# the last file written for each (marker set, replicate) and marks completion.

rule compare_IAblocks_with_fixedblocks:
    input:
        JULIA_ENV,
        "src/julia/comparison_block_vs_IA.jl",
        "src/julia/collection.jl",
        "data/rfmix_sequence.tsv",
        "data/rfmix_array.tsv",
        ref = "data/replicate_{r}/ancestral_{m}_10_10.vcf.gz",
        adm = "data/replicate_{r}/admixed_{m}_10_10.vcf.gz",
        anc = "data/replicate_{r}/global_ancestries.csv",
        tru = "data/replicate_{r}/local_ancestries.csv",
    output:
        "results/blocks/{m}_{r}_0.001_1.1.csv"
    threads: 10
    resources:
        mem_mb = lambda wc: _mem("localancestryjl", wc.m),
        runtime = lambda wc: 60 * _runtime("localancestryjl", wc.m),
    shell:
        "julia -t {threads} --project=. src/julia/comparison_block_vs_IA.jl {wildcards.m} {wildcards.r}"
