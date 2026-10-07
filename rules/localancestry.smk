
rule estimate_localancestryjl:
    input:
        "src/julia/setup.jl",
        "deps/localancestry.so",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries.csv"
    output:
        result = "results/localancestry/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
    threads: 24
    resources:
        mem_mb = 15000,
        runtime = 5
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
        "results/localancestry/rfmix/{r}_{m}_{g1}_{g2}.fb.tsv",
        "results/localancestry/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        "results/localancestry/rfmix/{r}_{m}_{g1}_{g2}.rfmix.Q",
        "results/localancestry/rfmix/{r}_{m}_{g1}_{g2}.sis.tsv"
    threads: 24
    resources:
        mem_mb = 100000,
        runtime = 720
    benchmark:
        "benchmarks/rfmix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        "rfmix/rfmix --n-threads={threads} -f {input.admixed} -r {input.ancestral} -m {input.globalancestries} -g {input.genomicmap} -o results/localancestry/rfmix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} --chromosome=chr22 "


rule estimate_flare:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        "results/localancestry/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = 20000,
        runtime = 30
    benchmark:
        "benchmarks/flare/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.0 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                out=results/localancestry/flare/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        """

rule estimate_flareclust:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        "results/localancestry/flareclust/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = 20000,
        runtime = 35
    benchmark:
        "benchmarks/flareclust/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.0 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                panel-probs=true \
                out=results/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        
            python3 flare/create_model_file.py \
            3 \
            results/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}.panels \
            results/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}

            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.0 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                model=results/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}.model \
                update-p=true \
                out=results/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        """

rule estimate_gnomix:
    input:
        "gnomix/config.yaml",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        ancestral_tbi = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz.tbi",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        admixed_tbi = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz.tbi",
        globalancestries = "data/replicate_{r}/global_ancestries_gnomix_{g1}.tsv",
        genomicmap = "data/rfmix_{m}.tsv",
    output:
        "results/localancestry/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp",
    threads: 24
    resources:
        mem_mb = 400000,
        runtime = 200
    benchmark:
        "benchmarks/gnomix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            OMP_NUM_THREADS=1
            OPENBLAS_NUM_THREADS=1

            PYTHONWARNINGS="ignore::FutureWarning" \
            python3 gnomix/gnomix.py \
            {input.admixed} \
            results/localancestry/gnomix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} \
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
        poplabels = "data/replicate_{r}/global_ancestries_rfmix.tsv",
    output:
        "results/localancestry/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt",
    threads: 24
    resources:
        mem_mb = 5000,
        runtime = 30
    benchmark:
        "benchmarks/recombmix/{r}_{m}_{g1}_{g2}.tsv"
    params:
        outprefix = "results/localancestry/recombmix/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
            Recomb-Mix/RecombMix -t {threads} -p {input.ancestral} -a {input.poplabels} -q {input.admixed} -g {input.genomicmap} -o results/localancestry/recombmix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        """

rule estimate_sparsepainter:
    input:
        reffile = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        targetfile = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        mapfile = "data/sparsepainter_{m}.tsv",
        popfile = "data/replicate_{r}/global_ancestries_rfmix.tsv",
        namefile = "data/replicate_{r}/individuals_admixed_gen{g2}_10.txt",
    output:
        prob = "results/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz",
    threads: 24
    resources:
        mem_mb = 400000,
        runtime = 180
    benchmark:
        "benchmarks/sparsepainter/{r}_{m}_{g1}_{g2}.tsv"
    params:
        out = "results/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
        
        # Expose LAPACK and BLAS
        export CMAKE_PREFIX_PATH="$CONDA_PREFIX${{CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}}"
        export LIBRARY_PATH="$CONDA_PREFIX/lib${{LIBRARY_PATH:+:$LIBRARY_PATH}}"
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"
        export CPATH="$CONDA_PREFIX/include${{CPATH:+:$CPATH}}"
        
        SparsePainter/SparsePainter \
            -reffile {input.reffile} \
            -targetfile {input.targetfile} \
            -mapfile {input.mapfile} \
            -popfile {input.popfile} \
            -namefile {input.namefile} \
            -out {params.out} \
            -prob \
            -probstore raw \
            -ncores {threads} 
        """

