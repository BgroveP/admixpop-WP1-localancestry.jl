
# Configfile
configfile: "config.yaml"

# Load variables from the config file
replicates = config["replicates"]   # r
generations = config["generations"] # g
populations = config["populations"] # p
markersets = config["markersets"]   # m

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
        "gnomix/gnomix.py"
    shell:
        r"""
        julia --project=. -e 'println("WOHOO")'
        """

# 
include: "rules/setup.smk"
