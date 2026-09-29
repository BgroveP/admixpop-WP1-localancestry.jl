
rule estimate_localancestryjl:
    input:
        "src/julia/setup.jl",
        "deps/localancestry.so",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries.csv"
    output:
        result = "results/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
    threads: 24
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.tsv"
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
        "rfmix/rfmix --n-threads={threads} -f {input.admixed} -r {input.ancestral} -m {input.globalancestries} -g {input.genomicmap} -o results/rfmix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} --chromosome=chr22"


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

rule estimate_gnomix:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        ancestral_tbi = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz.tbi",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        admixed_tbi = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz.tbi",
        globalancestries = "data/replicate_{r}/global_ancestries_gnomix_{}.ssv",
        genomicmap = "data/rfmix_{m}.tsv",
    output:
        "results/gnomix/{r}_{m}_{g1}_{g2}/chr22.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/gnomix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            python3 gnomix/gnomix.py \
            {input.admixed} \
            results/gnomix/{r}_{m}_{g1}_{g2} \
            chr22 \
            False \
            {input.genomicmap} \
            {input.ancestral} \
            {input.globalancestries} \
            gnomix/config.yaml
        """

        
rule estimate_recombmix:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        ancestral_tbi = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz.tbi",
        admixed_tbi = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz.tbi",
        genomicmap = "data/recombmix_{m}.tsv",
        poplabels = "data/replicate_{r}/global_ancestries_rfmix.txt",
    output:
        "results/recombmix/{r}_{m}_{g1}_{g2}/admix_inferred_ancestral_values_local.txt",
    threads: 1
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/recombmix/{r}_{m}_{g1}_{g2}.tsv"
    params:
        outprefix = "results/recombmix/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
            Recomb-Mix/RecombMix -t {threads} -p {input.ancestral} -a {input.poplabels} -q {input.admixed} -g {input.genomicmap} -o results/recombmix/{r}_{m}_{g1}_{g2}
        """

rule estimate_sparsepainter:
    input:
        ancestral_hap = "data/replicate_{r}/ancestral_{m}_{g1}_10.hap.gz",
        admixed_hap = "data/replicate_{r}/admixed_{m}_{g2}_10.hap.gz",
        reffile = "data/replicate_{r}/global_ancestries_sparsepainter.txt",
        tgtfile = "data/replicate_{r}/targets_{m}_{g2}_10.txt",
        inifile = "src/sparsepainter.ini",
    output:
        "results/sparsepainter/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = 30000,
        runtime = 180
    benchmark:
        "benchmarks/sparsepainter/{r}_{m}_{g1}_{g2}.tsv"
    params:
        outprefix = "results/sparsepainter/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
            SparsePainter \
                -h {input.ancestral_hap} \
                -a {input.admixed_hap} \
                -f {input.reffile} \
                -t {input.tgtfile} \
                -i {input.inifile} \
                -o {params.outprefix} \
                -c {threads} \
                > {output} 2> {log}
        """