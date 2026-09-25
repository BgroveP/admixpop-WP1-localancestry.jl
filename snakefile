import tarfile
import requests
from pathlib import Path

# Global variables
## Paths
DATA_DIR = "data"
ZENODO_FILE = f"{DATA_DIR}/raw/simulated_data.tgz"  
ZENODO_DESTINATION = f"{DATA_DIR}/extracted"

## URLs
ZENODO_URL = "https://zenodo.org/records/19544041/files/simulated_data.tgz?download=1"
FLARE_URL = "https://github.com/browning-lab/flare"
LOCALANCESTRY_URL = "https://github.com/BgroveP/LocalAncestry.jl"
RECOMBMIX_URL = "https://github.com/ucfcbb/Recomb-Mix"

## 
REPLICATES = range(1, 31)
GENERATIONS = [6, 8, 10]
POPULATIONS = ["ancestral", "admixed"]

# Download the compressed file from Zenodo
rule all:
    input:
        f"{ZENODO_DESTINATION}/replicate_30/admixed.vcf.gz",
        "Manifest.toml",
        "bin/flare.jar",
        "SparsePainter/SparsePainter",
        "Recomb-Mix/RecombMix"

rule zenodo_download:
    output:
        ZENODO_FILE  
    params:
        url = ZENODO_URL
    shell:
        "aria2c -x 15 -s 15 -o {output} {params.url}"

# Decompress the compressed file from Zenodo
rule zenodo_extract:
    output:
        expand(f"{ZENODO_DESTINATION}/replicate_{{i}}/{{name}}", i=REPLICATES, name=["admixed.vcf.gz", "admixed.vcf.gz.tbi", "ancestral.vcf.gz", "ancestral.vcf.gz.tbi", "global_ancestries.csv", "local_ancestries.csv"])
    input:
        ZENODO_FILE
    shell:
        "tar -xvzf {ZENODO_FILE} -C {ZENODO_DESTINATION}/"

# Setup , LocalAncestry.jl, FLARE, Recombmix, and RFMix
rule julia_setup:
    input:
        "src/setup.jl"
    output:
        "Manifest.toml",
        "Project.toml"
    shell:
        r"""
        julia --project=. -e 'using Pkg; Pkg.add(name="LocalAncestry", version="0.3.0")'
        """

rule flare_setup:
    output:
        "bin/flare.jar"
    shell:
        r"""
        # Clone the repository (remove existing directory if it exists)
        rm -rf flare
        git clone https://github.com/browning-lab/flare.git

        # Compile the Java file
        javac -cp flare/src/ -d flare/src/ flare/src/admix/AdmixMain.java

        # Create the bin directory and build the JAR
        mkdir -p bin && jar cfe bin/flare.jar admix/AdmixMain -C flare/src/ .
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

        #
        pixi add autoconf automake libtool

        #
        autoreconf --force --install # creates the configure script and all its dependencies
        ./configure                  # generates the Makefile
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
        pixi add openblas lapack
        rm -rf SparsePainter
        git clone https://github.com/YaolingYang/SparsePainter.git
        cd SparsePainter

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
        cd gnomix
        pip install -r requirements.txt
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
        "1kgenomes_all_chr22.vcf.gz".
        "1kgenomes_all_chr22.vcf.gz.tbi".
        "1kgenomes_all_individiuals.tsv"
    shell:
        r"""
        mkdir -p data/1kgenomes
        cd data/1kgenomes
        curl -o 1kgenomes_all_chr22.vcf.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
        curl -o 1kgenomes_all_chr22.vcf.gz.tbi https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
        curl -o 1kgenomes_all_individiuals.tsv https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        """

# Run ADMIXTURE
rule admixture:
    input:
        "data/1kgenomes/1kgenomes_all_chr22.vcf.gz".
        "data/1kgenomes/1kgenomes_all_chr22.vcf.gz.tbi".
        "data/1kgenomes/1kgenomes_all_individiuals.tsv",
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

# Create subsets
rule subset_vcfs:
    input:
        "data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005.Q"
    shell:
        r"""
            julia --project=. subset_vcfs.jl
        """

# Estimate local ancestries
rule estimate:
    input:
        "gnomix/gnomix",
        "rfmix/rfmix",
        "Manifest.toml",
        "SparsePainter/SparsePainter"
        "Recomb-Mix/RecombMix",
        expand(f"{ZENODO_DESTINATION}/replicate_{{i}}/{{name}}", i=REPLICATES, name=["admixed.vcf.gz", "admixed.vcf.gz.tbi", "ancestral.vcf.gz", "ancestral.vcf.gz.tbi", "global_ancestries.csv", "local_ancestries.csv"])
    output:
        ""
    shell:
        r"""
        
        """

# Collect output
