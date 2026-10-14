

def _runtime(tool, markerset):
    return config["estimatelimits"][tool][markerset]["runtime"]


def _mem(tool, markerset):
    return config["estimatelimits"][tool][markerset]["mem_mb"]


rule estimate_localancestryjl:
    input:
        "src/julia/setup.jl",
        "deps/localancestry.so",
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries.csv"
    output:
        result = config["paths"]["tmpdir"] + "/localancestry/localancestryjl/{r}_{m}_{g1}_{g2}_{t}.csv",
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
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_rfmix.tsv",
        genomicmap = "data/rfmix_{m}.tsv",
    output:
        config["paths"]["tmpdir"] + "/localancestry/rfmix/{r}_{m}_{g1}_{g2}.fb.tsv",
        config["paths"]["tmpdir"] + "/localancestry/rfmix/{r}_{m}_{g1}_{g2}.msp.tsv",
        config["paths"]["tmpdir"] + "/localancestry/rfmix/{r}_{m}_{g1}_{g2}.rfmix.Q",
        config["paths"]["tmpdir"] + "/localancestry/rfmix/{r}_{m}_{g1}_{g2}.sis.tsv"
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("rfmix", wc.m),
        runtime = lambda wc: _runtime("rfmix", wc.m),
    benchmark:
        "benchmarks/rfmix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        "rfmix/rfmix --n-threads={threads} -f {input.admixed} -r {input.ancestral} -m {input.globalancestries} -g {input.genomicmap} -o {resources.tmpdir}/localancestry/rfmix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} --chromosome=chr22 "


rule estimate_flare:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        config["paths"]["tmpdir"] + "/localancestry/flare/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("flare", wc.m),
        runtime = lambda wc: _runtime("flare", wc.m),
    benchmark:
        "benchmarks/flare/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.0 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                out={resources.tmpdir}/localancestry/flare/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        """


rule estimate_flareclust:
    input:
        ancestral = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        admixed = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        globalancestries = "data/replicate_{r}/global_ancestries_flare.ssv",
        genomicmap = "data/flare_{m}.ssv",
    output:
        config["paths"]["tmpdir"] + "/localancestry/flareclust/{r}_{m}_{g1}_{g2}.anc.vcf.gz",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("flareclust", wc.m),
        runtime = lambda wc: _runtime("flareclust", wc.m),
    benchmark:
        "benchmarks/flareclust/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.0 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                panel-probs=true \
                out={resources.tmpdir}/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        
            python3 flare/create_model_file.py \
            3 \
            {resources.tmpdir}/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}.panels \
            {resources.tmpdir}/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}

            java -jar flare/flare.jar \
                min-mac=0 min-maf=0.0 \
                ref={input.ancestral} ref-panel={input.globalancestries} \
                gt={input.admixed} map={input.genomicmap} \
                model={resources.tmpdir}/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2}.model \
                update-p=true \
                out={resources.tmpdir}/localancestry/flareclust/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
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
        config["paths"]["tmpdir"] + "/localancestry/gnomix/{r}_{m}_{g1}_{g2}/query_results.msp",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("gnomix", wc.m),
        runtime = lambda wc: _runtime("gnomix", wc.m),
    benchmark:
        "benchmarks/gnomix/{r}_{m}_{g1}_{g2}.tsv"
    shell:
        r"""
            OMP_NUM_THREADS=1
            OPENBLAS_NUM_THREADS=1

            PYTHONWARNINGS="ignore::FutureWarning" \
            python3 gnomix/gnomix.py \
            {input.admixed} \
            {resources.tmpdir}/localancestry/gnomix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} \
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
        config["paths"]["tmpdir"] + "/localancestry/recombmix/{r}_{m}_{g1}_{g2}/inferred_local_ancestral_values.txt",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("recombmix", wc.m),
        runtime = lambda wc: _runtime("recombmix", wc.m),
    benchmark:
        "benchmarks/recombmix/{r}_{m}_{g1}_{g2}.tsv"
    params:
        outprefix = config["paths"]["tmpdir"] + "/localancestry/recombmix/{r}_{m}_{g1}_{g2}",
    shell:
        r"""
            Recomb-Mix/RecombMix -t {threads} -p {input.ancestral} -a {input.poplabels} -q {input.admixed} -g {input.genomicmap} -o {resources.tmpdir}/localancestry/recombmix/{wildcards.r}_{wildcards.m}_{wildcards.g1}_{wildcards.g2} 
        """


rule estimate_sparsepainter:
    input:
        reffile = "data/replicate_{r}/ancestral_{m}_{g1}_10.vcf.gz",
        targetfile = "data/replicate_{r}/admixed_{m}_{g2}_10.vcf.gz",
        mapfile = "data/sparsepainter_{m}.tsv",
        popfile = "data/replicate_{r}/global_ancestries_sparsepainter_{g1}.tsv",
        namefile = "data/replicate_{r}/individuals_admixed_gen{g2}_10.txt",
    output:
        prob = config["paths"]["tmpdir"] + "/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}_prob.txt.gz",
    threads: 24
    resources:
        mem_mb = lambda wc: _mem("sparsepainter", wc.m),
        runtime = lambda wc: _runtime("sparsepainter", wc.m),
    benchmark:
        "benchmarks/sparsepainter/{r}_{m}_{g1}_{g2}.tsv"
    params:
        out = config["paths"]["tmpdir"] + "/localancestry/sparsepainter/{r}_{m}_{g1}_{g2}",
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