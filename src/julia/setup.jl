
# Install the Julia packages used by the pipeline and build a system image for LocalAncestry.jl.
# Usage: julia --project=. src/julia/setup.jl <LocalAncestry version> <QGIO commit>
# (called by the Snakemake rule `setup_julia`; versions are set in config.yaml)
# Run it through a cluster profile so it is compiled on a compute node, e.g.
#   pixi run snakemake --profile profiles/ghpcaudk setup_julia

using Pkg

localancestry_version = length(ARGS) >= 1 ? ARGS[1] : "0.3.1"
qgio_commit = length(ARGS) >= 2 ? ARGS[2] : "c196f98fbffc372913de72b4c22af4207cba2f0f"

# Do not precompile after each Pkg.add; everything is precompiled once below
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"

# Install all packages in a single resolve: registered packages, LocalAncestry.jl at the
# pinned version and a VCF parser from GitHub
Pkg.add([
    [PackageSpec(name = p) for p in ["DataFrames", "CSV", "StatsBase", "PackageCompiler", "CodecZlib",
                                     "Combinatorics", "AlgebraOfGraphics", "CairoMakie"]];
    PackageSpec(name = "LocalAncestry", version = localancestry_version);
    PackageSpec(url = "https://github.com/BgroveP/QGIO.jl", rev = qgio_commit)
])

# Make sure everything in the manifest is installed, then precompile all packages once,
# here, for the current JULIA_CPU_TARGET. Otherwise the first pipeline jobs would all
# precompile the same packages in parallel (slow, and fragile on a shared depot).
println("JULIA_CPU_TARGET = ", get(ENV, "JULIA_CPU_TARGET", "(not set: caches are only valid for this CPU type)"))
println("Parallel precompile tasks = ", get(ENV, "JULIA_NUM_PRECOMPILE_TASKS", "(default)"))
Pkg.instantiate()
Pkg.precompile(strict = true)   # fail the job if any package does not precompile

# Compile system image
using PackageCompiler
mkpath("deps")
create_sysimage(
    [:CSV, :DataFrames, :LocalAncestry, :CodecZlib];
    sysimage_path = "deps/localancestry.so",
    # Portable across x86-64 CPU generations: the image may be built on a different
    # node type than the one that later uses it
    cpu_target = PackageCompiler.default_app_cpu_target(),
    precompile_execution_file = "deps/precompile_workload.jl"
)
