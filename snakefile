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
        "plots/accuracy_ancgen.png",
        "plots/accuracy_software.png",
        "plots/sensitivities.png",
        "plots/accuracy_blocksize.png"
    resources:
        mem_mb = 100000,
        runtime = 180
    output:
        "finished.txt"
    threads: 24
    shell:
        r"""
        touch finished.txt
        """

# 
include: "rules/1. dataprep.smk"
include: "rules/2. softwaresetup.smk"
include: "rules/3. localancestry.smk"
include: "rules/4. benchmark.smk"
include: "rules/5. iablocksVSfixedblocks.smk"
include: "rules/6. plots.smk"
