
# Software setup. Versions/commits are pinned in config.yaml (section `software`).
sw = config["software"]

# Julia packages (incl. LocalAncestry.jl) and a precompiled system image
rule setup_julia:
    input:
        "src/julia/setup.jl",
        "deps/precompile_workload.jl"
    output:
        "deps/localancestry.so"
    params:
        version = sw["localancestry"],
        qgio = sw["qgio_commit"]
    # Package precompilation and system image code generation run in parallel.
    # Set the task counts explicitly: Julia otherwise sizes them from the node's total
    # core count, not from the cores SLURM allocated to the job.
    threads: 16
    resources:
        mem_mb = 32000,
        runtime = 180
    shell:
        """
        export JULIA_NUM_PRECOMPILE_TASKS={threads}
        export JULIA_IMAGE_THREADS={threads}
        julia --project=. -t {threads} src/julia/setup.jl {params.version} {params.qgio}
        """

rule flare_setup:
    output:
        "flare/flare.jar"
    params:
        commit = sw["flare_commit"]
    shell:
        r"""
        # Clone the repository (remove existing directory if it exists)
        rm -rf flare
        git clone https://github.com/browning-lab/flare.git
        cd flare
        git checkout {params.commit}

        # Compile the Java file and build the JAR
        javac -cp src/ -d src/ src/admix/AdmixMain.java
        jar cfe flare.jar admix/AdmixMain -C src/ .
        """

rule rfmix_setup:
    output:
        "rfmix/rfmix"
    params:
        commit = sw["rfmix_commit"]
    shell:
        r"""
        rm -rf rfmix
        git clone https://github.com/slowkoni/rfmix.git
        cd rfmix
        git checkout {params.commit}

        # Prepare the make file and compile
        autoreconf --force --install
        ./configure
        make
        """

rule recombmix_setup:
    output:
        "Recomb-Mix/RecombMix"
    params:
        commit = sw["recombmix_commit"]
    shell:
        r"""
        rm -rf Recomb-Mix
        git clone https://github.com/ucfcbb/Recomb-Mix.git
        git -C Recomb-Mix checkout {params.commit}

        # Build in a separate pixi environment with its own compiler and boost
        cd Recomb-Mix
        pixi init .
        pixi add gcc_linux-64 gxx_linux-64 make boost-cpp
        pixi run bash -c '
            sed -i "s|^CXXFLAGS = .*|CXXFLAGS = -I\$\$CONDA_PREFIX/include -std=c++17 -Os -fopenmp|" Makefile
            sed -i "s|^LIBS = .*|LIBS = -L\$\$CONDA_PREFIX/lib -Wl,-rpath,\$\$CONDA_PREFIX/lib -lboost_iostreams -fopenmp|" Makefile
            make
        '
        """

rule sparsepainter_setup:
    output:
        "SparsePainter/SparsePainter"
    params:
        commit = sw["sparsepainter_commit"]
    shell:
        r"""
        rm -rf SparsePainter
        git clone https://github.com/YaolingYang/SparsePainter.git
        cd SparsePainter
        git checkout {params.commit}

        export CMAKE_PREFIX_PATH="$CONDA_PREFIX${{CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}}"
        export LIBRARY_PATH="$CONDA_PREFIX/lib${{LIBRARY_PATH:+:$LIBRARY_PATH}}"
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"
        export CPATH="$CONDA_PREFIX/include${{CPATH:+:$CPATH}}"

        make
        """

rule gnomix_setup:
    output:
        "gnomix/gnomix.py",
        "gnomix/config.yaml"
    params:
        commit = sw["gnomix_commit"]
    shell:
        r"""
        rm -rf gnomix
        git clone https://github.com/AI-sandbox/gnomix
        git -C gnomix checkout {params.commit}
        pip install -r gnomix/requirements.txt

        # Choose the "large" model
        sed -i 's/^\(\s*inference:\s*\)\(#\)/\1"large" \2/' gnomix/config.yaml

        # Bugfix for outdated function argument when using "large" model
        sed -i 's/padding_mode="reflection"/padding_mode="reflect"/' gnomix/src/Smooth/cnn.py
        """

rule admixture_setup:
    output:
        "admixture/admixture"
    shell:
        r"""
        curl -fsSLO https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz
        tar xf admixture_linux-1.4.0.tar.gz
        rm -f admixture_linux-1.4.0.tar.gz
        mkdir -p admixture
        mv -u admixture_linux-1.4.0/* ./admixture/
        rm -rf admixture_linux-1.4.0
        """
