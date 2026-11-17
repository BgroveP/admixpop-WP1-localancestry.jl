
# Download the simulated data from Zenodo
rule zenodo_download:
    output:
        config["paths"]["zenodopath"]
    params:
        url = config["paths"]["zenodourl"]
    threads: 15
    shell:
        r"""
        aria2c -x {threads} -s {threads} -o {output} {params.url}
        """

# Unzip the simulated data
rule zenodo_extract:
    input:
        config["paths"]["zenodopath"]
    output:
        expand("data/replicate_{r}/{file}", r=replicates, file=config["paths"]["zenodocontent"])
    shell:
        r"""
        mkdir -p data
        tar -xzf {input} -C data/
        """

# Create individual files for subsetting vcfs across generations
rule individual_files:
    input:
        JULIA_ENV,
        expand("data/replicate_{r}/global_ancestries.csv", r=replicates),
        "src/julia/create_individual_files.jl"
    output:
        expand("data/replicate_{r}/individuals_{c}_gen{g}_10.txt", r=replicates, g=generations, c=cohorts),
        expand("data/replicate_{r}/global_ancestries_rfmix.tsv", r=replicates),
        expand("data/replicate_{r}/global_ancestries_gnomix_{g}.tsv", r=replicates, g=generations),
        expand("data/replicate_{r}/global_ancestries_flare.ssv", r=replicates)
    params:
        generations = ",".join(str(g) for g in generations)
    shell:
        "julia --project=. src/julia/create_individual_files.jl {params.generations}"

# SparsePainter needs a population file restricted to the reference individuals.
# The individuals are the same for every marker set, so the first/smallest configured one is used.
rule individual_files_sparsepainter:
    input:
        JULIA_ENV,
        poplabels = "data/replicate_{r}/global_ancestries.csv",
        jlscript = "src/julia/write_sparsepainter_popfile.jl",
        vcfpath = "data/replicate_{r}/ancestral_" + markersets[0] + "_{g1}_10.vcf.gz"
    output:
        outpath = "data/replicate_{r}/global_ancestries_sparsepainter_{g1}.tsv"
    shell:
        "julia --project=. {input.jlscript} {input.vcfpath} {input.poplabels} {output.outpath}"

# Download HapMap3 marker positions (used to define the "array" marker set) and lift them over to hg38
rule hapmap3_download_and_liftover:
    output:
        "data/hapmap3.hg38.bed"
    shell:
        r"""
        rm -rf data/tmp/hapmap3
        mkdir -p data/tmp/hapmap3
        cd data/tmp/hapmap3

        # download liftOver and chain files
        curl -fsSLO http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
        curl -fsSLO http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
        curl -fsSLO http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        chmod 0744 liftOver

        # download hapmap3 data in plink format
        curl -fsSLO https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
        bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2

        # map -> bed, then hg18 -> hg19 -> hg38
        awk -F'\t' -v OFS='\t' '{{ print "chr"$1, $4, $4 }}' hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map > hapmap3.hg18.bed
        ./liftOver hapmap3.hg18.bed hg18ToHg19.over.chain.gz hapmap3.hg19.bed unmapped_hg18_to_hg19.log
        ./liftOver hapmap3.hg19.bed hg19ToHg38.over.chain.gz ../../hapmap3.hg38.bed unmapped_hg19_to_hg38.log
        """

# Download 1000 Genomes chromosome 22
rule x1Kgenomes_download:
    output:
        "data/1kgenomes_all_chr22.vcf.gz",
        "data/1kgenomes_all_chr22.vcf.gz.tbi",
        "data/1kgenomes_all_individuals.tsv"
    params:
        base = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
    shell:
        r"""
        curl -fL -o {output[0]} {params.base}
        curl -fL -o {output[1]} {params.base}.tbi
        curl -fL -o {output[2]} https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        """

# Keep 1000 Genomes AFR/EAS/EUR individuals with >= 99% ADMIXTURE ancestry in their own super-population
rule admixture:
    input:
        JULIA_ENV,
        "data/1kgenomes_all_chr22.vcf.gz",
        "data/1kgenomes_all_chr22.vcf.gz.tbi",
        "data/1kgenomes_all_individuals.tsv",
        "src/julia/filter_1kgenomes.jl",
        "src/julia/filter_1kgenomes2.jl",
        "admixture/admixture"
    output:
        "data/all_individuals.tsv",
        "data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz",
        "data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz.tbi"
    shell:
        r"""
        # ADMIXTURE needs a large stack; not all systems allow raising it
        ulimit -s unlimited || echo "Warning: could not set unlimited stack size" >&2

        # Download IGSR individuals
        curl -f \
            -X POST \
            -H 'Content-Type: application/x-www-form-urlencoded' \
            --data-urlencode 'json={{"fields":["name","populations.superpopulationCode"],"column_names":["individual","superpopulation"]}}' \
            'https://www.internationalgenome.org/api/beta/sample/_search/igsr_samples.tsv' \
            -o data/all_individuals.tsv

        # Subset for AFR, EAS, and EUR, and find intersect between individuals files
        julia --project=. src/julia/filter_1kgenomes.jl

        # Subset the VCF and index
        bcftools view -S data/1kgenomes_AFR_EAS_EUR_individuals.txt -Oz -o data/1kgenomes_AFR_EAS_EUR_chr22.vcf.gz data/1kgenomes_all_chr22.vcf.gz
        tabix -p vcf data/1kgenomes_AFR_EAS_EUR_chr22.vcf.gz

        # Omit loci with no variation
        bcftools view -i 'MAF>0.005' data/1kgenomes_AFR_EAS_EUR_chr22.vcf.gz -Oz -o data/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz
        tabix -p vcf data/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz

        # Run admixture
        cd data/
        plink --vcf 1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz --make-bed --out 1kgenomes_AFR_EAS_EUR_chr22_maf005
        ../admixture/admixture 1kgenomes_AFR_EAS_EUR_chr22_maf005.bed 3
        cd ..
        julia --project=. src/julia/filter_1kgenomes2.jl
        bcftools view -S data/1kgenomes_AFR_EAS_EUR_individuals_ga99p.txt -Oz -o data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz data/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz
        tabix -p vcf data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz
        """

# Genomic maps in each tool's format, for both marker sets
rule genomic_maps:
    input:
        JULIA_ENV,
        "data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz",
        "data/hapmap3.hg38.bed",
        "src/julia/write_genomic_maps.jl"
    output:
        expand("data/flare_{ms}.ssv", ms=ALL_MARKERSETS),
        expand("data/rfmix_{ms}.tsv", ms=ALL_MARKERSETS),
        expand("data/recombmix_{ms}.tsv", ms=ALL_MARKERSETS),
        expand("data/sparsepainter_{ms}.tsv", ms=ALL_MARKERSETS),
        expand("data/bcftools_{ms}.tsv", ms=ALL_MARKERSETS)
    shell:
        r"""
        julia --project=. src/julia/write_genomic_maps.jl
        """

# Subset the simulated VCFs to a marker set and a range of generations
rule subset_vcfs:
    input:
        vcf = "data/replicate_{r}/{c}.vcf.gz",
        vcftbi = "data/replicate_{r}/{c}.vcf.gz.tbi",
        ind = "data/replicate_{r}/individuals_{c}_gen{g}_10.txt",
        gmap = "data/bcftools_{m}.tsv"
    output:
        vcf = "data/replicate_{r}/{c}_{m}_{g}_10.vcf.gz",
        vcftbi = "data/replicate_{r}/{c}_{m}_{g}_10.vcf.gz.tbi"
    threads: 4
    resources:
        runtime = 180
    shell:
        r"""
        bcftools view --threads {threads} -T {input.gmap} -S {input.ind} {input.vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """
