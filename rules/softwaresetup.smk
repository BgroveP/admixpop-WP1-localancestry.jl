
# Setup , LocalAncestry.jl, FLARE, Recombmix, and RFMix
rule julia_setup:
    input:
        srcfile = "src/julia/setup.jl"
    output:
        manifestfile = "Manifest.toml",
        projectfile = "Project.toml"
    shell:
        r"""
        rm -f {output.manifestfile}
        rm -f {output.projectfile}
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
        pixi init .
        pixi add gcc_linux-64 gxx_linux-64 make boost-cpp
        pixi install
        pixi shell -c '
        # Compile program
                sed -i 's|^CXXFLAGS = .*|CXXFLAGS = -I$(CONDA_PREFIX)/include -std=c++17 -Os -fopenmp|' Makefile
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

