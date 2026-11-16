
# Local ancestry inference with each tool.
# Raw output goes to LA_DIR (= <tmpdir>/localancestry; tmpdir is set in the Snakemake profile).
# Runtime/memory come from `estimatelimits` in config.yaml.

rule estimate_localancestryjl:
    input:
        "deps/localancestry.so",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries.csv"
    output:
        result = LA_DIR + "/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("localancestryjl", wc.m),
        runtime = lambda wc: _runtime("localancestryjl", wc.m),
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
                threshold={wildcards.t},
                maf=0.000
            )

            CSV.write("{output.result}", lajl)
        '
        """


rule estimate_rfmix:
    input:
        "rfmix/rfmix",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_rfmix.tsv",
        genomicmap = "data/rfmix_{m}.tsv",
    output:
        LA_DIR + "/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        LA_DIR + "/rfmix/{r}_{m}_{g1}_{g2}.rfmix.Q",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("rfmix", wc.m),
        runtime = lambda wc: _runtime("rfmix", wc.m),
    benchmark:
        "benchmarks/rfmix/{r}_{m}_{g1}_{g2}.tsv"
    params:
        out = LA_DIR + "/rfmix/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
        rfmix/rfmix --n-threads={threads} -f {input.admixed} -r {input.ancestral} -m {input.globalancestries} -g {input.genomicmap} -o {params.out} --chromosome=chr22
        # Large files that are not used downstream
        rm -f {params.out}.fb.tsv {params.out}.sis.tsv
        """


rule estimate_flare:
    input:
        "flare/flare.jar",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        LA_DIR + "/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("flare", wc.m),
        runtime = lambda wc: _runtime("flare", wc.m),
    params:
        heap_mb = lambda wc: int(_mem("flare", wc.m) * 0.75),
        out = LA_DIR + "/flare/{r}_{m}_{g1}_{g2}",
    benchmark:
        "benchmarks/flare/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
        java -Xmx{params.heap_mb}m -jar flare/flare.jar \
            min-mac=0 min-maf=0.0 \
            ref={input.ancestral} ref-panel={input.globalancestries} \
            gt={input.admixed} map={input.genomicmap} \
            out={params.out}
        """


rule estimate_flareclust:
    input:
        "flare/flare.jar",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        LA_DIR + "/flareclust/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("flareclust", wc.m),
        runtime = lambda wc: _runtime("flareclust", wc.m),
    params:
        heap_mb = lambda wc: int(_mem("flareclust", wc.m) * 0.75),
        out = LA_DIR + "/flareclust/{r}_{m}_{g1}_{g2}",
    benchmark:
        "benchmarks/flareclust/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
        java -Xmx{params.heap_mb}m -jar flare/flare.jar \
            min-mac=0 min-maf=0.0 \
            ref={input.ancestral} ref-panel={input.globalancestries} \
            gt={input.admixed} map={input.genomicmap} \
            panel-probs=true \
            out={params.out}

        python3 flare/create_model_file.py 3 {params.out}.panels {params.out}

        java -Xmx{params.heap_mb}m -jar flare/flare.jar \
            min-mac=0 min-maf=0.0 \
            ref={input.ancestral} ref-panel={input.globalancestries} \
            gt={input.admixed} map={input.genomicmap} \
            model={params.out}.model \
            update-p=true \
            out={params.out}
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
        LA_DIR + "/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("gnomix", wc.m),
        runtime = lambda wc: _runtime("gnomix", wc.m),
    params:
        out = LA_DIR + "/gnomix/{r}_{m}_{g1}_{g2}",
    benchmark:
        "benchmarks/gnomix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        PYTHONWARNINGS="ignore::FutureWarning" \
        python3 gnomix/gnomix.py \
            {input.admixed} \
            {params.out} \
            chr22 \
            False \
            {input.genomicmap} \
            {input.ancestral} \
            {input.globalancestries} \
            gnomix/config.yaml

        # Large intermediate files that are not used downstream
        rm -rf {params.out}/generated_data
        """


rule estimate_recombmix:
    input:
        "Recomb-Mix/RecombMix",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        ancestral_tbi = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz.tbi",
        admixed_tbi = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz.tbi",
        genomicmap = "data/recombmix_{m}.tsv",
        poplabels = "data/replicate_{r}/global_ancestries_rfmix.tsv",
    output:
        LA_DIR + "/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("recombmix", wc.m),
        runtime = lambda wc: _runtime("recombmix", wc.m),
    benchmark:
        "benchmarks/recombmix/{r}_{m}_{g1}_{g2}.tsv"
    params:
        out = LA_DIR + "/recombmix/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
        Recomb-Mix/RecombMix -t {threads} -p {input.ancestral} -a {input.poplabels} -q {input.admixed} -g {input.genomicmap} -o {params.out}
        """


rule estimate_sparsepainter:
    input:
        "SparsePainter/SparsePainter",
        reffile = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        targetfile = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        mapfile = "data/sparsepainter_{m}.tsv",
        popfile = "data/replicate_{r}/global_ancestries_sparsepainter_{g1}.tsv",
        namefile = "data/replicate_{r}/individuals_admixed_gen{g2}_10.txt",
    output:
        prob = LA_DIR + "/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("sparsepainter", wc.m),
        runtime = lambda wc: _runtime("sparsepainter", wc.m),
    benchmark:
        "benchmarks/sparsepainter/{r}_{m}_{g1}_{g2}.tsv"
    params:
        out = LA_DIR + "/sparsepainter/{r}_{m}_{g1}_{g2}",
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
