import re

# Configfile
configfile: "config.yaml"

# Load variables from the config file
replicates = config["replicates"]       # r
generations = config["generations"]     # g
populations = config["populations"]     # p
cohorts = config["cohorts"]     # p
markersets = config["markersets"]       # m
lathresholds = config["lathresholds"]  # t

def wc_constraint(values):
    alts = sorted((re.escape(str(v)) for v in values), key=len, reverse=True)
    return "|".join(alts)

wildcard_constraints:
    r = wc_constraint(replicates),
    g = wc_constraint(generations),
    g1 = wc_constraint(generations),
    g2 = wc_constraint(generations),
    p = wc_constraint(populations),
    c = wc_constraint(cohorts),
    m = wc_constraint(markersets),
    t = wc_constraint(lathresholds)


#
zenodocontent = config["paths"]["zenodocontent"]

# Download the compressed file from Zenodo
rule comparison:
    input:
        "Project.toml",
        expand("results/accuracy/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations, t = lathresholds),
        expand("results/accuracy/rfmix/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/accuracy/flare/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),
        #expand("results/accuracy/flareclust/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/accuracy/gnomix/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/accuracy/recombmix/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/accuracy/sparsepainter/{r}_{m}_{g1}_{g2}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations)
    resources:
        mem_mb = 100000,
        runtime = 180
    output:
        "results/summary/RAM_use_mean.csv",
        "results/summary/RAM_use_max.csv",
        "results/summary/time_use_mean.csv",
        "results/summary/time_use_max.csv",
        "plots/Accuracy_chromosomewide.png",
        "plots/Accuracy_chromosomewide_localancestryjl_thresholds.png"
    threads: 24
    shell:
        r"""
        julia --project=. src/julia/result_benchmark.jl
        julia --project=. src/julia/result_accuracies.jl
        """

# 
include: "rules/softwaresetup.smk"
include: "rules/dataprep.smk"
include: "rules/localancestry.smk"
include: "rules/benchmark.smk"
