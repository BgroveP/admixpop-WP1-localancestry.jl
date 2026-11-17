
# LocalAncestry.jl accuracy across a grid of maximum block sizes and IA thresholds.
# Usage: julia -t <threads> --project=. src/julia/comparison_block_vs_IA.jl <markerset> <replicate>
# Writes results/blocks/<markerset>_<replicate>_<blocksize>_<threshold>.csv
# (the last file written, *_0.001_1.1.csv, is the Snakemake completion marker).

using CSV
using DataFrames
using LocalAncestry
using Base.Threads
using StatsBase

include("collection.jl")

# Input
m = ARGS[1]
r = ARGS[2]

# Grid. `maxblockfraction` must be in (0, 1], so a tiny offset is added below.
blocksizes = 0.0:0.0005:0.003
thresholds = 0.3:0.1:1.1

# Paths
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"
refpath = "data/replicate_$(r)/ancestral_$(m)_10_10.vcf.gz"
admpath = "data/replicate_$(r)/admixed_$(m)_10_10.vcf.gz"
ancpath = "data/replicate_$(r)/global_ancestries.csv"
trupath = "data/replicate_$(r)/local_ancestries.csv"
outfile(b, t) = "results/blocks/$(m)_$(r)_$(b)_$(t).csv"
mkpath("results/blocks")

# Genomic map with array indicator
gmap = CSV.read(gmappath, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2 = CSV.read(gmappath2, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2.array .= true
leftjoin!(gmap, gmap2, on = names(gmap))
gmap.array[ismissing.(gmap.array)] .= false

# True local ancestries
tru = CSV.read(trupath, DataFrame)
tru.basepairs = UnitRange{Int}.(parse.(Int, replace.(tru.basepairs, r"\:.+" => "")), parse.(Int, replace.(tru.basepairs, r".+\:" => "")))

function run_one(b, t)
    ofile = outfile(b, t)
    isfile(ofile) && return nothing   # resume: skip finished combinations
    thisgmap = copy(gmap)
    est = localancestry(refpath, admpath, ancpath, threshold = t, maf = 0.0, maxblockfraction = b + 10^-9)
    accuracies = evaluate(est, tru, thisgmap)
    thisgmap.accuracy = accuracies.pcorrect
    thisgmap.incorrect = accuracies.pincorrect
    thisgmap.unassigned = accuracies.punassigned
    println("$(m) replicate $(r), block size $(b), threshold $(t): ", mean(thisgmap.accuracy[Bool.(thisgmap.array)]))

    # Only array loci are used in the comparison
    deleteat!(thisgmap, .~thisgmap.array)
    CSV.write(ofile, thisgmap)
    return nothing
end

# Run all combinations except the completion marker, then the marker last
marker = (0.001, 1.1)
combos = [(b, t) for b in blocksizes for t in thresholds if !(b ≈ marker[1] && t ≈ marker[2])]
@threads for i in eachindex(combos)
    run_one(combos[i]...)
end
run_one(marker...)
