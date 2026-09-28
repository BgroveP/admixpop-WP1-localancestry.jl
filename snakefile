
# Configfile
configfile: "config.yaml"

# Load variables from the config file
replicates = config["replicates"]       # r
generations = config["generations"]     # g
populations = config["populations"]     # p
markersets = config["markersets"]       # m
lathresholds = config["lathresholds"]   # t

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
        expand("data/replicate_{r}/{p}_{m}_{g}_10.vcf.gz", r = replicates, m = markersets, p =populations, g = generations),
        expand("results/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.tsv", r = replicates, m = markersets, p =populations, g1 = generations, g2 = generations, t = lathresholds),
        expand("results/rfmix/{r}_{m}_{g1}_{g2}.fb.tsv", r = replicates, m = markersets, p =populations, g1 = generations, g2 = generations),
        expand("results/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz", r = replicates, m = markersets, p =populations, g1 = generations, g2 = generations)
    shell:
        r"""
        julia --project=. -e 'println("WOHOO")'
        """

# 
include: "rules/softwaresetup.smk"
include: "rules/dataprep.smk"
include: "rules/localancestry.smk"
