
# Per-locus accuracy of each tool against the simulated (true) local ancestries,
# plus summaries of accuracy and computational resources.

# Inputs shared by every accuracy script
ACCURACY_COMMON = [
    JULIA_ENV,
    "src/julia/collection.jl",
    "data/rfmix_sequence.tsv",   # genomic maps define the evaluated loci
    "data/rfmix_array.tsv",
]

rule accuracy_flare:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_flare.jl",
        "src/julia/read_flare.jl",
        est = LA_DIR + "/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/flare/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        "julia --project=. src/julia/accuracy_flare.jl {input.est} {input.tru} {output.path}"

rule accuracy_flareclust:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_flareclust.jl",
        "src/julia/read_flare.jl",
        est = LA_DIR + "/flareclust/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/flareclust/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        "julia --project=. src/julia/accuracy_flareclust.jl {input.est} {input.tru} {output.path}"

rule accuracy_rfmix:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_rfmix.jl",
        "src/julia/read_rfmix.jl",
        est = LA_DIR + "/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        tru = "data/replicate_{r}/local_ancestries.csv",
        loci = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",  # loci used to resolve window ends
    output:
        path = "results/accuracy/rfmix/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        "julia --project=. src/julia/accuracy_rfmix.jl {input.est} {input.tru} {output.path} {input.loci}"

rule accuracy_recombmix:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_recombmix.jl",
        "src/julia/read_recombmix.jl",
        est = LA_DIR + "/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/recombmix/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        "julia --project=. src/julia/accuracy_recombmix.jl {input.est} {input.tru} {output.path}"

rule accuracy_localancestryjl:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_localancestryjl.jl",
        est = LA_DIR + "/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv"
    resources:
        mem_mb = 10000
    shell:
        "julia --project=. src/julia/accuracy_localancestryjl.jl {input.est} {input.tru} {output.path}"

rule accuracy_gnomix:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_gnomix.jl",
        "src/julia/read_rfmix.jl",
        est = LA_DIR + "/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp",
        tru = "data/replicate_{r}/local_ancestries.csv",
        loci = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",  # loci used to resolve window ends
    output:
        path = "results/accuracy/gnomix/{r}_{m}_{g1}_{g2}.csv"
    resources:
        mem_mb = 10000
    shell:
        "julia --project=. src/julia/accuracy_gnomix.jl {input.est} {input.tru} {output.path} {input.loci}"

rule accuracy_sparsepainter:
    input:
        ACCURACY_COMMON,
        "src/julia/accuracy_sparsepainter.jl",
        "src/julia/read_sparsepainter.jl",
        est = LA_DIR + "/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz",
        tru = "data/replicate_{r}/local_ancestries.csv"
    output:
        path = "results/accuracy/sparsepainter/{r}_{m}_{g1}_{g2}.csv"
    threads: 10
    resources:
        mem_mb = 40000
    shell:
        "julia -t {threads} --project=. src/julia/accuracy_sparsepainter.jl {input.est} {input.tru} {output.path}"


# --- Summaries ---
# Arguments to the summary scripts: <comma-separated tools> <default LocalAncestry.jl threshold>
SUMMARY_ARGS = ",".join(tools) + " " + str(DEFAULT_THRESHOLD)

rule benchmark_accuracy:
    input:
        JULIA_ENV,
        "src/julia/result_accuracies.jl",
        tool_files("results/accuracy/{tool}/{r}_{m}_{g1}_{g2}{suffix}.csv")
    output:
        "results/summary/accuracy.csv",
        "results/summary/accuracy_sds.csv",
        "results/summary/missing_accuracy.csv"
    params:
        args = SUMMARY_ARGS
    threads: 10
    resources:
        mem_mb = 50000
    shell:
        "julia -t {threads} --project=. src/julia/result_accuracies.jl {params.args}"

rule benchmark_computation_resources:
    input:
        JULIA_ENV,
        "src/julia/result_benchmark.jl",
        tool_files("benchmarks/{tool}/{r}_{m}_{g1}_{g2}{suffix}.tsv")
    output:
        "results/summary/RAM_use_mean.csv",
        "results/summary/time_use_mean.csv",
        "results/summary/RAM_use_sds.csv",
        "results/summary/time_use_sds.csv",
        "results/summary/resources_table_latex.csv"
    params:
        args = SUMMARY_ARGS
    shell:
        "julia --project=. src/julia/result_benchmark.jl {params.args}"
