
rule plot_accuracy:
    input:
        "results/summary/accuracy.csv",
        "results/summary/accuracy_sds.csv",
    output:
        "plots/accuracy_blocksize.png",
    shell:
        r"""
            julia --project=. src/julia/result_plot.jl
        """

rule plot_blocksize:
    input:
        "results/summary/accuracy.csv",
        "results/summary/accuracy_sds.csv",
    output:
        "plots/accuracy_ancgen.png",
        "plots/accuracy_software.png"
    shell:
        r"""
            julia --project=. src/julia/result_block.jl
        """

rule plot_sensitivity:
    input:
        "results/summary/accuracy.csv",
        "results/summary/RAM_use_mean.csv",
        "results/summary/time_use_mean.csv",
    output:
        "plots/sensitivities.png"
    shell:
        r"""
            julia --project=. src/julia/result_sensitivity.jl
        """
