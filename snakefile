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
rule all:
    input:
        "Project.toml",
        "flare/flare.jar",
        "rfmix/rfmix",
        "Recomb-Mix/RecombMix",
        "SparsePainter/SparsePainter",
        "gnomix/gnomix.py",
        expand("data/replicate_{r}/{c}_{m}_{g}_10.vcf.gz", r = replicates, m = markersets, c = cohorts, g = generations),
        expand("results/localancestry/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv", r = replicates, m = markersets, g1 = generations, g2 = generations, t = lathresholds),
        expand("results/localancestry/rfmix/{r}_{m}_{g1}_{g2}.fb.tsv", r = replicates, m = ["array"], g1 = generations, g2 = generations),
        expand("results/localancestry/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/localancestry/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/localancestry/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt", r = replicates, m = markersets, g1 = generations, g2 = generations),
        expand("results/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz", r = replicates, m = markersets, g1 = generations, g2 = generations)
    shell:
        r"""
        julia --project=. -e 'println("WOHOO")'
        """

# 
include: "rules/softwaresetup.smk"
include: "rules/dataprep.smk"
include: "rules/localancestry.smk"
