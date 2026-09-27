
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
        "checkpoint.txt"

# 
include: "rules/setup.smk"

# Download the simulated data from Zenodo
rule zenodo_download:
    output:
        config["paths"]["zenodopath"]
    params:
        url = config["paths"]["zenodourl"]
    threads: config["threads"]["high"]
    shell:
        r"""
        mkdir -p $(dirname {output}
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

# Setup , LocalAncestry.jl, FLARE, Recombmix, and RFMix
rule julia_setup:
    input:
        srcfile = "src/julia/setup.jl"
    output:
        "Manifest.toml",
        "Project.toml"
    shell:
        r"""
        julia --project=. {input.srcfile}
        """

rule flare_setup:
    output:
        "flare/flare.jar"
    shell:
        r"""
        # Clone the repository (remove existing directory if it exists)
        rm -rf flare
        git clone https://github.com/browning-lab/flare.git
        cd flare
        git checkout d6a6c8fd4479ebdb68529824acd11d2b6bbac4a1

        # Compile the Java file
        javac -cp src/ -d src/ src/admix/AdmixMain.java

        # Create the bin directory and build the JAR
        jar cfe flare.jar admix/AdmixMain -C src/ .
        """

rule rfmix_setup:
    output:
        "rfmix/rfmix"
    shell:
        r"""
        # Clone the repository (remove existing directory if it exists)
        rm -rf rfmix
        git clone https://github.com/slowkoni/rfmix.git
        cd rfmix
        git checkout 09599c19c69f81571c7a7aab4cd2d847a9e0633a

        # Prepare the make file and compile
        autoreconf --force --install 
        ./configure                  
        make
        """

rule recombmix_setup:
    output:
        "Recomb-Mix/RecombMix"
    shell:
        r"""
        # Clone the repository (remove existing directory if it exists)
        rm -rf Recomb-Mix
        git clone https://github.com/ucfcbb/Recomb-Mix.git
        cd Recomb-Mix
        git checkout 41616840ccdfd456c7645c4ee651450a890994d5 

        # Initiate the pixi environment for recombmix
        pixi init recombmix
        pixi add gcc_linux-64 gxx_linux-64 make boost-cpp
        pixi install
        pixi shell -c '

        # Compile program
        sed -i 's|^CXXFLAGS = |CXXFLAGS = -Ic/include/boost |' Makefile
        sed -i 's|^LIBS = .*|LIBS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib -lboost_iostreams -fopenmp|' Makefile
        make
        '
        """

rule sparsepainter_setup:
    output:
        "SparsePainter/SparsePainter"
    shell:
        r"""
        rm -rf SparsePainter
        git clone https://github.com/YaolingYang/SparsePainter.git
        cd SparsePainter
        git checkout ed91149c5ecf795daba28711e08138f8abb00f68

        export CMAKE_PREFIX_PATH="$CONDA_PREFIX${{CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}}"
        export LIBRARY_PATH="$CONDA_PREFIX/lib${{LIBRARY_PATH:+:$LIBRARY_PATH}}"
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"
        export CPATH="$CONDA_PREFIX/include${{CPATH:+:$CPATH}}"

        make
        """

rule gnomix_setup:
    output:
        "gnomix/gnomix.py"
    shell:
        r"""
        rm -rf gnomix
        git clone https://github.com/AI-sandbox/gnomix
        pip install -r gnomix/requirements.txt
        """

# Create individual files for subsetting vcfs across generations
rule individual_files:
    input:
        expand(f"{ZENODO_DESTINATION}/replicate_{{i}}/global_ancestries.csv", i=REPLICATES),
        "src/create_individual_files.jl"
    output:
        expand(f"{ZENODO_DESTINATION}/replicate_{{i}}/individuals_{{p}}_gen{{g}}_10.txt", i=REPLICATES, g = GENERATIONS, p = POPULATIONS)
    shell:
        "julia --project=. src/create_individual_files.jl"

# Download hapmap3 map
rule hapmap3_download_and_liftover:
    output:
        "data/maps/hapmap3.hg38.bed"
    shell:
        r"""
        mkdir -p tmp data/maps
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
        awk -F'\t' -v OFS='\t' '{ print "chr"$1, $4, $4 }' hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map > hapmap3.hg18.bed
        ./liftOver hapmap3.hg18.bed hg18ToHg19.over.chain.gz hapmap3.hg19.bed unmapped_hg18_to_hg19.log
        ./liftOver hapmap3.hg19.bed hg19ToHg38.over.chain.gz ../data/maps/hapmap3.hg38.bed unmapped_hg19_to_hg38.log
        """

# Download 1K genomes 
rule x1Kgenomes_download:
    output:
        "1kgenomes_all_chr22.vcf.gz",
        "1kgenomes_all_chr22.vcf.gz.tbi",
        "1kgenomes_all_individuals.tsv"
    shell:
        r"""
        mkdir -p data/1kgenomes
        cd data/1kgenomes
        curl -o 1kgenomes_all_chr22.vcf.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
        curl -o 1kgenomes_all_chr22.vcf.gz.tbi https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
        curl -o 1kgenomes_all_individuals.tsv https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        """

# Run ADMIXTURE
rule admixture:
    input:
        "data/1kgenomes/1kgenomes_all_chr22.vcf.gz",
        "data/1kgenomes/1kgenomes_all_chr22.vcf.gz.tbi",
        "data/1kgenomes/1kgenomes_all_individuals.tsv",
        "src/filter_1kgenomes.jl",
        "src/filter_1kgenomes2.jl"
    output:
        "data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz",
        "data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz.tbi"
    shell:
        r"""
        # Download IGSR individuals
        curl \
          -X POST \
          -H 'Content-Type: application/x-www-form-urlencoded' \
          --data-urlencode 'json={"fields":["name","populations.superpopulationCode"],"column_names":["individual","superpopulation"]}' \
          'https://www.internationalgenome.org/api/beta/sample/_search/igsr_samples.tsv' \
          -o data/1kgenomes/all_individuals.tsv

        # Subset for AFR, EAS, and EUR, and find intersect between individuals files
        julia --project=. src/filter_1kgenomes.jl

        # Subset the VCF and index
        bcftools view -S data/1kgenomes/1kgenomes_AFR_EAS_EUR_individuals.txt -Oz -o data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22.vcf.gz data/1kgenomes/1kgenomes_all_chr22.vcf.gz
        tabix -p vcf data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22.vcf.gz

        # Omit loci with no variation
        bcftools view -i 'MAF>0.005' data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22.vcf.gz -Oz -o data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz
        tabix -p vcf data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz

        # Run admixture
        cd data/1kgenomes/
        plink --vcf 1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz --make-bed --out 1kgenomes_AFR_EAS_EUR_chr22_maf005
        admixture 1kgenomes_AFR_EAS_EUR_chr22_maf005.bed 3 
        cd ../..
        julia --project=. src/filter_1kgenomes2.jl
        bcftools view -S data/1kgenomes/1kgenomes_AFR_EAS_EUR_individuals_ga99p.txt -Oz -o data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005.vcf.gz
        tabix -p vcf data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz
        """

rule genomic_maps:
    input:
        "data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz",
        "data/maps/hapmap3.hg38.bed",
        "src/write_genomic_maps.jl"
    output:
        expand(f"data/maps/flare_{{m}}.ssv", m = MARKERSETS),
        expand(f"data/maps/rfmix_{{m}}.tsv", m = MARKERSETS),
        expand(f"data/maps/recombmix_{{m}}.tsv", m = MARKERSETS),
        expand(f"data/maps/bcftools_{{m}}.bed", m = MARKERSETS)
    shell:
        r"""
        julia --project=. src/write_genomic_maps.jl 
        """

rule subset_vcfs:
    input:
        vcf = "data/extracted/replicate_{r}/{p}.vcf.gz",
        vcftbi = "data/extracted/replicate_{r}/{p}.vcf.gz.tbi",
        ind = "data/extracted/replicate_{r}/individuals_{p}_gen{g}.txt",
        gmap = "data/maps/bcftools_{m}.bed"
    output:
        vcf = "data/extracted/replicate_{r}/{p}_{m}_{g}.vcf.gz",
        vcftbi = "data/extracted/replicate_{r}/{p}_{m}_{g}.vcf.gz.tbi"
    shell:
        r"""
            bcftools view -R {input.gmap} -S {input.ind} {input.vcf} -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        """
    
rule estimate_localancestryjl:
    input:
        ancestral = f"{ZENODO_DESTINATION}/replicate_{{i}}/ancestral_{{markerset}}_{{ancgens}}_10.vcf.gz",
        admixed = f"{ZENODO_DESTINATION}/replicate_{{i}}/admixed_{{markerset}}_{{admgens}}_10.vcf.gz",
        globalancestries = f"{ZENODO_DESTINATION}/replicate_{{i}}/localancestryjl_globalancestry.csv",
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
