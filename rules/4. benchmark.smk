
rule accuracy_flare:
    input:
        "src/julia/accuracy_flare.jl",
        "src/julia/read_flare.jl",
        "src/julia/collection.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/flare/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        r"""
            julia --project=. src/julia/accuracy_flare.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_flareclust:
    input:
        "src/julia/accuracy_flare.jl",
        "src/julia/read_flare.jl",
        "src/julia/collection.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/flareclust/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/flareclust/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        r"""
            julia --project=. src/julia/accuracy_flareclust.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_rfmix:
    input:
        "src/julia/accuracy_rfmix.jl",
        "src/julia/read_rfmix.jl",
        "src/julia/collection.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/rfmix/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        r"""
            julia --project=. src/julia/accuracy_rfmix.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_recombmix:
    input:
        "src/julia/accuracy_recombmix.jl",
        "src/julia/read_recombmix.jl",
        "src/julia/collection.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/recombmix/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        r"""
            julia --project=. src/julia/accuracy_recombmix.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_localancestryljl:
    input:
        "src/julia/accuracy_localancestryjl.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv"
    resources:
        mem_mb = 10000
    shell:
        r"""
            julia --project=. src/julia/accuracy_localancestryjl.jl {input.est} {input.tru} {output.path}
        """


rule accuracy_gnomix:
    input:
        "src/julia/accuracy_gnomix.jl",
        "src/julia/read_rfmix.jl",
        "src/julia/collection.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/gnomix/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        r"""
            julia --project=. src/julia/accuracy_gnomix.jl {input.est} {input.tru} {output.path}
        """

rule accuracy_sparsepainter:
    input:
        "src/julia/accuracy_sparsepainter.jl",
        "src/julia/collection.jl",
        est = config["paths"]["tmpdir"] + "/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/sparsepainter/{r}_{m}_{g1}_{g2}.csv"
    threads:
        10
    resources:
        mem_mb = 40000
    shell:
        r"""
            julia -t {threads} --project=. src/julia/accuracy_sparsepainter.jl {input.est} {input.tru} {output.path}
        """
rule benchmark_accuracy:
    input:
        expand("results/accuracy/gnomix/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("results/accuracy/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations, t = lathresholds),    
        expand("results/accuracy/rfmix/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("results/accuracy/recombmix/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("results/accuracy/sparsepainter/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("results/accuracy/flare/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations)    
    output:
        "results/summary/accuracy.csv",
        "results/summary/accuracy_sds.csv"
    threads:
        10
    shell:
        r"""
            julia -t {threads} --project=. src/julia/result_accuracies.jl
        """

rule benchmark_computation_resources:
    input:
        expand("benchmarks/gnomix/{r}_{m}_{g1}_{g2}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("benchmarks/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations, t = lathresholds),    
        expand("benchmarks/rfmix/{r}_{m}_{g1}_{g2}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("benchmarks/recombmix/{r}_{m}_{g1}_{g2}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("benchmarks/sparsepainter/{r}_{m}_{g1}_{g2}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations),    
        expand("benchmarks/flare/{r}_{m}_{g1}_{g2}.tsv", r = replicates, m = markersets, g1 = generations, g2 = generations)    
    output:
        "results/summary/RAM_use_mean.csv",
        "results/summary/time_use_mean.csv"
    shell:
        r"""
            julia --project=. src/julia/result_benchmark.jl
        """
