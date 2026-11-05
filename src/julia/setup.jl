
# Load the Pkg manager
using Pkg

# Install official packages
Pkg.add(["DataFrames", "CSV", "StatsBase", "PackageCompiler", "CodecZlib", "Combinatorics", "AlgebraOfGraphics", "CairoMakie"])
Pkg.add(name="LocalAncestry", version="0.4.0")

# Install a VCF parser from github
Pkg.add(url= "https://github.com/BgroveP/QGIO.jl", rev="c196f98fbffc372913de72b4c22af4207cba2f0f")

# Compile system image
using PackageCompiler
create_sysimage(
    [:CSV, :DataFrames, :LocalAncestry, :CodecZlib];
    sysimage_path = "deps/localancestry.so",
    precompile_execution_file = "deps/precompile_workload.jl"
)

