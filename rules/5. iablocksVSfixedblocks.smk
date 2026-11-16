
rule compare_IAblocks_with_fixedblocks:
    input:
        expand("results/accuracy/localancestryjl/{r}_{m}_{g1}_{g2}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations, t = lathresholds),    
    output:
        expand("results/blocks/{m}_{r}_0.001_1.1.csv", r = replicates, m = markersets)
    threads: 10
    shell:
        r"""
            julia -t {threads} --project=. src/julia/comparison_block_vs_IA.jl
        """
