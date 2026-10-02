
rule accuracy_flare:
    input:
        "src/julia/accuracy_flare.jl",
        "src/julia/read_flare.jl",
        "src/julia/collection.jl",
        est = "results/localancestry/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/flare/{r}_{m}_{g1}_{g2}.csv"
    shell:
        r"""
            julia --project=. src/julia/accuracy_flare.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_rfmix:
    input:
        "src/julia/accuracy_rfmix.jl",
        "src/julia/read_rfmix.jl",
        "src/julia/collection.jl",
        est = "results/localancestry/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/rfmix/{r}_{m}_{g1}_{g2}.csv"
    shell:
        r"""
            julia --project=. src/julia/accuracy_rfmix.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_recombmix:
    input:
        "src/julia/accuracy_recombmix.jl",
        "src/julia/read_recombmix.jl",
        "src/julia/collection.jl",
        est = "results/localancestry/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/recombmix/{r}_{m}_{g1}_{g2}.csv"
    shell:
        r"""
            julia --project=. src/julia/accuracy_recombmix.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_localancestryljl:
    input:
        "src/julia/accuracy_localancestryjl.jl",
        est = "results/localancestry/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv"
    shell:
        r"""
            julia --project=. src/julia/accuracy_localancestryjl.jl {input.est} {input.tru} {output.path}
        """


rule accuracy_gnomix:
    input:
        "src/julia/accuracy_gnomix.jl",
        "src/julia/read_rfmix.jl",
        "src/julia/collection.jl",
        est = "results/localancestry/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/gnomix/{r}_{m}_{g1}_{g2}.csv"
    shell:
        r"""
            julia --project=. src/julia/accuracy_gnomix.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_sparsepainter:
    input:
        "src/julia/accuracy_sparsepainter.jl",
        "src/julia/collection.jl",
        est = "results/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/sparsepainter/{r}_{m}_{g1}_{g2}.csv"
    shell:
        r"""
            julia --project=. src/julia/accuracy_sparsepainter.jl {input.est} {input.tru} {output.path}
        """
