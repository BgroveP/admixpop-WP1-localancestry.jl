
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
        expand("data/replicate_{r}/{file}", r=replicates, file=zenodocontent)
    shell:
        r"""
        tar -xzf {input} -C data/
        """

# Create individual files for subsetting vcfs across generations
rule individual_files:
    input:
        expand("data/replicate_{r}/global_ancestries.csv", r=replicates),
        "src/julia/create_individual_files.jl"
    output:
        expand("data/replicate_{r}/individuals_{c}_gen{g}_10.txt", r=replicates, g = generations, c = cohorts),
        expand("data/replicate_{r}/global_ancestries_rfmix.tsv", r=replicates),
        expand("data/replicate_{r}/global_ancestries_gnomix_{g}.tsv", r=replicates, g = generations),
        expand("data/replicate_{r}/global_ancestries_flare.ssv", r=replicates)
    shell:
        "julia --project=. src/julia/create_individual_files.jl"

# Download hapmap3 map
rule hapmap3_download_and_liftover:
    output:
        "data/hapmap3.hg38.bed"
    shell:
        r"""
        mkdir -p tmp data
        cd tmp

        # download database and script
        wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        chmod 0744 liftOver 

        # download hapmap3 data in plink format
        curl -O https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
        bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2

        # Homemade
        awk -F'\t' -v OFS='\t' '{{ print "chr"$1, $4, $4 }}' hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map > hapmap3.hg18.bed
        ./liftOver hapmap3.hg18.bed hg18ToHg19.over.chain.gz hapmap3.hg19.bed unmapped_hg18_to_hg19.log
        ./liftOver hapmap3.hg19.bed hg19ToHg38.over.chain.gz ../data/hapmap3.hg38.bed unmapped_hg19_to_hg38.log
        """

# Download 1K genomes 
rule x1Kgenomes_download:
    output:
        "data/1kgenomes_all_chr22.vcf.gz",
        "data/1kgenomes_all_chr22.vcf.gz.tbi",
        "data/1kgenomes_all_individuals.tsv"
    shell:
        r"""
        cd data
        curl -o 1kgenomes_all_chr22.vcf.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
        curl -o 1kgenomes_all_chr22.vcf.gz.tbi https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
        curl -o 1kgenomes_all_individuals.tsv https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        """

# Run ADMIXTURE
rule admixture:
    input:
        "data/1kgenomes_all_chr22.vcf.gz",
        "data/1kgenomes_all_chr22.vcf.gz.tbi",
        "data/1kgenomes_all_individuals.tsv",
        "src/julia/filter_1kgenomes.jl",
        "src/julia/filter_1kgenomes2.jl",
        "Manifest.toml",
        "Project.toml"
    output:
        "data/all_individuals.tsv",
        "data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz",
        "data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz.tbi"
    shell:
        r"""
        # Download IGSR individuals
        curl \
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
        admixture 1kgenomes_AFR_EAS_EUR_chr22_maf005.bed 3 
        cd ..
        julia --project=. src/julia/filter_1kgenomes2.jl
        bcftools view -S data/1kgenomes_AFR_EAS_EUR_individuals_ga99p.txt -Oz -o data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz data/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz
        tabix -p vcf data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz
        """

rule genomic_maps:
    input:
        "data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz",
        "data/hapmap3.hg38.bed",
        "src/julia/write_genomic_maps.jl"
    output:
        expand("data/flare_{m}.ssv", m = markersets),
        expand("data/rfmix_{m}.tsv", m = markersets),
        expand("data/recombmix_{m}.tsv", m = markersets),
        expand("data/sparsepainter_{m}.tsv", m = markersets),
        expand("data/bcftools_{m}.tsv", m = markersets)
    shell:
        r"""
        julia --project=. src/julia/write_genomic_maps.jl 
        """

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
