
rule estimate_localancestryjl:
    input:
        ancestral = "data/replicate_{i}/ancestral_{markerset}_{ancgens}_10.vcf.gz",
        admixed = "data/replicate_{i}/admixed_{markerset}_{admgens}_10.vcf.gz",
        globalancestries = "data/replicate_{i}/localancestryjl_globalancestry.csv",
    output:
        result = "results/localancestry/localancestryjl_{i}_{markerset}_{ancgens}_{admgens}_{lathreshold}.tsv",
    threads: 24
    resources:
        mem_mb = 30000,
        slurm_account = "admixpop"
    benchmark:
        "results/localancestryjl_{i}_{markerset}_{ancgens}_{admgens}_{lathreshold}.tsv"
    shell:
        r"""
        mkdir -p results

        julia --project=. -t auto -e '
            using CSV 
            using DataFrames
            using LocalAncestry

            lajl = localancestry(
                "{input.ancestral}",
                "{input.admixed}",
                "{input.globalancestries}",
                threshold={wildcards.lathreshold}
            )

            CSV.write("{output.result}", lajl)
        '
        """

# Estimate local ancestries with LocalAncestry.jl
rule estimate:
    input:
        "gnomix/gnomix.py",
        "rfmix/rfmix",
        "Manifest.toml",
        "SparsePainter/SparsePainter",
        "Recomb-Mix/RecombMix",
        expand(f"{ZENODO_DESTINATION}/replicate_{{i}}/{{p}}_{{m}}_{{g}}_10.vcf.gz", i=REPLICATES, p=["ancestral", "admixed"], m=["array", "sequence"], g = ["6", "8", "10"]),
        expand(f"data/maps/flare_{{m}}.ssv", m=["array", "sequence"]),
        expand(f"data/maps/rfmix_{{m}}.tsv", m=["array", "sequence"]),
        expand(f"data/maps/recombmix_{{m}}.tsv", m=["array", "sequence"])
    output:
        "checkpoint.txt"
    shell:
        r"""
        touch checkpoint.txt
        """

# Collect output
