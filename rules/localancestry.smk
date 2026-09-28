
rule estimate_localancestryjl:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries.csv",
    output:
        result = "results/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
    threads: 24
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv"
    shell:
        r"""
        julia -J deps/localancestry.so --startup-file=no --project=. -t {threads} -e '
            using CSV 
            using DataFrames
            using LocalAncestry

            lajl = localancestry(
                "{input.ancestral}",
                "{input.admixed}",
                "{input.globalancestries}",
                threshold={wildcards.t}
            )

            CSV.write("{output.result}", lajl)
        '
        """

rule estimate_rfmix:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_rfmix.tsv",
        genomicmap = "data/rfmix_{m}.tsv",
    output:
        "results/rfmix/{r}_{m}_{g1}_{g2}.fb.tsv",
        "results/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        "results/rfmix/{r}_{m}_{g1}_{g2}.rfmix.Q",
        "results/rfmix/{r}_{m}_{g1}_{g2}.sis.tsv"
    threads: 24
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/rfmix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        "rfmix/rfmix --n-threads={threads} -f {input.admixed} -r {input.ancestral} -m {input.globalancestries} -g {input.genomicmap} -o results/rfmix_{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} --chromosome=chr22"


rule estimate_flare:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        "results/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/flare/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.001 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                out=results/flare/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}
        """

