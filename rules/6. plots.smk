
# Figures. Each rule's output matches the files written by its script.

rule plot_accuracy:
    input:
        JULIA_ENV,
        "src/julia/result_plot.jl",
        "results/summary/accuracy.csv",
        "results/summary/accuracy_sds.csv",
    output:
        "plots/accuracy_ancgen.png",
        "plots/accuracy_software.png",
    params:
        args = SUMMARY_ARGS
    shell:
        "julia --project=. src/julia/result_plot.jl {params.args}"

rule plot_blocksize:
    input:
        JULIA_ENV,
        "src/julia/result_block.jl",
        expand("results/blocks/{m}_{r}_0.001_1.1.csv", m=markersets, r=replicates),
    output:
        "plots/accuracy_blocksize.png",
        "plots/accuracy_blocksize.pdf",
    shell:
        "julia --project=. src/julia/result_block.jl"

rule plot_sensitivity:
    input:
        JULIA_ENV,
        "src/julia/result_sensitivity.jl",
        "results/summary/accuracy.csv",
        "results/summary/RAM_use_mean.csv",
        "results/summary/time_use_mean.csv",
    output:
        "plots/sensitivities.png",
        "plots/sensitivities.pdf",
    params:
        args = SUMMARY_ARGS
    shell:
        "julia --project=. src/julia/result_sensitivity.jl {params.args}"
