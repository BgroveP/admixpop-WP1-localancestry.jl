
# Real data: LocalAncestry.jl on an admixed 1000 Genomes Phase 3 population (default ACB),
# genome-wide, compared with the published 1000 Genomes Phase 3 local ancestry calls
# (RFMix v1.5.4; technical/working/20140818_ancestry_deconvolution, hg19).
#
# Genotypes are the Phase 3 release VCFs (GRCh37), i.e. the same build as the published
# calls, so no liftover is needed. Per chromosome:
#   biallelic SNPs, MAF >= realdata.maf  ->  reference / target VCFs  ->  LocalAncestry.jl
#   ->  concordance with the published calls per individual.
# The genetic length of each chromosome (LocalAncestry.jl `totcm`) is taken from the cM
# columns of the published bed files (HapMap genetic map used by 1000 Genomes).

RDCFG = config["realdata"]
RD_TARGET = RDCFG["target"]
RD_MODE = RDCFG["deconvolution_mode"]
RD_CHROMS = RDCFG["chromosomes"]
RD = "data/realdata"
RD_BEDS = f"{RD}/deconvolution/{RD_TARGET}_{RD_MODE}_bed_files"
RD_LA = LA_DIR + "/realdata/localancestryjl"

wildcard_constraints:
    chrom = wc_constraint(RD_CHROMS)


rule realdata_panel:
    output:
        f"{RD}/panel.tsv"
    params:
        url = RDCFG["panelurl"]
    shell:
        "curl -fL --retry 3 -o {output} {params.url}"


# Published local ancestry calls: one bed file per haplotype (<sample>_A.bed, <sample>_B.bed)
rule realdata_deconvolution:
    output:
        directory(RD_BEDS)
    params:
        url = RDCFG["deconvolutionurl"].format(target=RD_TARGET),
        zip = f"{RD}/deconvolution/{RD_TARGET}_phase3_ancestry_deconvolution.zip",
        extracted = f"{RD}/deconvolution/{RD_TARGET}_extracted",
        mode = RD_MODE,
    shell:
        r"""
        curl -fL --retry 3 -o {params.zip} {params.url}
        rm -rf {params.extracted}
        python3 -m zipfile -e {params.zip} {params.extracted}

        # Expected layout: <POP>/<mode>/bed_files/ (see the release README)
        beddir=$(find {params.extracted} -type d -path "*/{params.mode}/bed_files" | head -n 1)
        if [ -z "$beddir" ]; then
            echo "Could not find {params.mode}/bed_files in {params.zip}. Contents:" >&2
            find {params.extracted} -maxdepth 3 >&2
            exit 1
        fi
        rm -rf {output}
        mv "$beddir" {output}
        rm -rf {params.extracted}
        echo "Published calls: $(ls {output} | wc -l) bed files"
        """


# Reference individuals (by population) and target individuals with both genotypes and published calls
rule realdata_samples:
    input:
        panel = f"{RD}/panel.tsv",
        beds = RD_BEDS,
    output:
        target = f"{RD}/samples_target.txt",
        reference = f"{RD}/samples_reference.txt",
        ancestries = f"{RD}/ancestries.csv",
        all = f"{RD}/samples_all.txt",
    params:
        target = RD_TARGET,
        popmap = ",".join(f"{pop}:{anc}" for anc, pops in RDCFG["references"].items() for pop in pops),
    shell:
        r"""
        awk -F'\t' -v t={params.target} 'NR > 1 && $2 == t {{ print $1 }}' {input.panel} | sort > {output.target}.panel
        ls {input.beds} | sed -n 's/_A\.bed$//p' | sort > {output.target}.calls
        comm -12 {output.target}.panel {output.target}.calls > {output.target}
        rm {output.target}.panel {output.target}.calls
        n=$(wc -l < {output.target})
        echo "{params.target} individuals with genotypes and published calls: $n"
        if [ "$n" -eq 0 ]; then echo "No {params.target} individuals in both the panel and the published calls" >&2; exit 1; fi

        awk -F'\t' -v map={params.popmap} '
            BEGIN {{ n = split(map, a, ","); for (i = 1; i <= n; i++) {{ split(a[i], b, ":"); m[b[1]] = b[2] }}; print "individual,population" }}
            NR > 1 && ($2 in m) {{ print $1 "," m[$2] }}' {input.panel} > {output.ancestries}
        tail -n +2 {output.ancestries} | cut -d, -f1 > {output.reference}
        echo "Reference individuals: $(wc -l < {output.reference})"
        cat {output.reference} {output.target} > {output.all}
        """


# Genetic length (cM) of each chromosome, from the published bed files
rule realdata_chromosome_lengths:
    input:
        RD_BEDS
    output:
        f"{RD}/chromosome_lengths_cm.tsv"
    shell:
        r"""
        cat {input}/*.bed | awk -F'\t' '
            {{ c = $1; sub(/^chr/, "", c)
               if (!(c in lo) || $5 < lo[c]) lo[c] = $5
               if (!(c in hi) || $6 > hi[c]) hi[c] = $6 }}
            END {{ for (c in hi) print c "\t" hi[c] - lo[c] }}' | sort -k1,1n > {output}
        """


# Download one chromosome, keep biallelic SNPs above the MAF threshold, split into reference and target
rule realdata_vcf:
    input:
        all = f"{RD}/samples_all.txt",
        reference = f"{RD}/samples_reference.txt",
        target = f"{RD}/samples_target.txt",
    output:
        reference = f"{RD}/vcf/reference_chr{{chrom}}.vcf.gz",
        reference_tbi = f"{RD}/vcf/reference_chr{{chrom}}.vcf.gz.tbi",
        target = f"{RD}/vcf/target_chr{{chrom}}.vcf.gz",
        target_tbi = f"{RD}/vcf/target_chr{{chrom}}.vcf.gz.tbi",
    params:
        url = RDCFG["vcfurl"],   # {chrom} is filled in by Snakemake
        raw = f"{RD}/vcf/raw_chr{{chrom}}.vcf.gz",
        filtered = f"{RD}/vcf/filtered_chr{{chrom}}.vcf.gz",
        maf = RDCFG["maf"],
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 240,
    shell:
        r"""
        curl -fL --retry 3 -o {params.raw} {params.url}

        bcftools view --threads {threads} -S {input.all} -m2 -M2 -v snps {params.raw} -Ou \
            | bcftools norm -d all -Ou \
            | bcftools view --threads {threads} -i 'MAF>{params.maf}' -Oz -o {params.filtered}
        rm -f {params.raw}

        bcftools view --threads {threads} -S {input.reference} {params.filtered} -Oz -o {output.reference}
        bcftools view --threads {threads} -S {input.target} {params.filtered} -Oz -o {output.target}
        tabix -p vcf {output.reference}
        tabix -p vcf {output.target}
        rm -f {params.filtered}
        """


rule realdata_localancestryjl:
    input:
        "deps/localancestry.so",
        reference = f"{RD}/vcf/reference_chr{{chrom}}.vcf.gz",
        target = f"{RD}/vcf/target_chr{{chrom}}.vcf.gz",
        ancestries = f"{RD}/ancestries.csv",
        lengths = f"{RD}/chromosome_lengths_cm.tsv",
    output:
        RD_LA + f"/{RD_TARGET}_chr{{chrom}}.csv"
    params:
        threshold = DEFAULT_THRESHOLD,
    threads: 24
    resources:
        mem_mb = RDCFG["limits"]["mem_mb"],
        runtime = RDCFG["limits"]["runtime"],
    benchmark:
        f"benchmarks/realdata/localancestryjl/{RD_TARGET}_chr{{chrom}}.tsv"
    shell:
        r"""
        julia -J deps/localancestry.so --startup-file=no --project=. -t {threads} -e '
            using CSV
            using DataFrames
            using LocalAncestry

            lengths = CSV.read("{input.lengths}", DataFrame, header = ["chromosome", "cm"], types = [String, Float64])
            totcm = only(lengths.cm[lengths.chromosome .== "{wildcards.chrom}"])
            println("Chromosome {wildcards.chrom}: genetic length ", totcm, " cM")

            la = localancestry(
                "{input.reference}",
                "{input.target}",
                "{input.ancestries}",
                chromosome = "{wildcards.chrom}",
                threshold = {params.threshold},
                maf = 0.0,
                totcm = totcm
            )

            CSV.write("{output}", la)
        '
        """


# Per-individual agreement with the published calls on one chromosome
rule realdata_concordance:
    input:
        JULIA_ENV,
        "src/julia/realdata_concordance.jl",
        est = RD_LA + f"/{RD_TARGET}_chr{{chrom}}.csv",
        target = f"{RD}/vcf/target_chr{{chrom}}.vcf.gz",
        beds = RD_BEDS,
        samples = f"{RD}/samples_target.txt",
    output:
        f"results/realdata/{RD_TARGET}_chr{{chrom}}_individuals.csv"
    resources:
        mem_mb = 16000,
    shell:
        "julia --project=. src/julia/realdata_concordance.jl {input.est} {input.target} {input.beds} {input.samples} {wildcards.chrom} {output}"


# Genome-wide summary of concordance and resource use
rule realdata_summary:
    input:
        JULIA_ENV,
        "src/julia/realdata_summary.jl",
        expand("results/realdata/{target}_chr{chrom}_individuals.csv", target=RD_TARGET, chrom=RD_CHROMS),
        expand("benchmarks/realdata/localancestryjl/{target}_chr{chrom}.tsv", target=RD_TARGET, chrom=RD_CHROMS),
    output:
        "results/realdata/summary.csv",
        "results/realdata/individuals.csv",
        "results/realdata/resources.csv",
    params:
        target = RD_TARGET,
        chroms = ",".join(str(c) for c in RD_CHROMS),
    shell:
        "julia --project=. src/julia/realdata_summary.jl {params.target} {params.chroms}"
