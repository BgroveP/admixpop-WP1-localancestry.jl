
# Load the library
using Pkg
Pkg.update()
Pkg.add(["CSV", "DataFrames", "StatsBase"])
Pkg.add(url = "https://github.com/BgroveP/LocalAncestry.jl")

using CSV
using DataFrames
using LocalAncestry
using Base.threads
using StatsBase

# Input
blocksizes = 0.0:0.0005:0.003
thresholds = reverse(0.2:0.05:1.1)
replicates = rand(1:30, 30)
markersets = reverse(["array", "sequence"])

# Constants
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"
include("collection.jl")

# Set up the genomic map
gmap = CSV.read(gmappath, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2 = CSV.read(gmappath2, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2.array .= true
leftjoin!(gmap, gmap2, on=names(gmap))
gmap.array[ismissing.(gmap.array)] .= false

# Get N
todos = ["results/blocks/$(m)_$(r)_$(b)_$(t).csv" for m in markersets, r in replicates, b in blocksizes, t in thresholds if !isfile("results/blocks/$(m)_$(r)_$(b)_$(t).csv")]

#
for m in markersets
    for r in replicates
    
    # Set paths
    refpath = "data/replicate_$(r)/ancestral_$(m)_10_10.vcf.gz"
    admpath = "data/replicate_$(r)/admixed_$(m)_10_10.vcf.gz"
    ancpath = "data/replicate_$(r)/global_ancestries.csv"
    trupath = "data/replicate_$(r)/local_ancestries.csv"

    # Read true local ancestries
    tru = CSV.read(trupath, DataFrame)
    tru.basepairs = UnitRange{Int}.(parse.(Int, replace.(tru.basepairs, r"\:.+" => "")), parse.(Int, replace.(tru.basepairs, r".+\:" => "")))

    @threads for b in blocksizes
        for t in thresholds
        ofile = "results/blocks/$(m)_$(r)_$(b)_$(t).csv"
            if !isfile(ofile)
            thisgmap = copy(gmap)
            est = localancestry(refpath, admpath, ancpath, threshold = t, maf = 0.0, maxblockfraction = b + 10^-9)
            accuracies = evaluate(est, tru, thisgmap)
            thisgmap.accuracy = accuracies.pcorrect
            thisgmap.incorrect = accuracies.pincorrect
            thisgmap.unassigned = accuracies.punassigned
            println(mean(thisgmap.accuracy[Bool.(thisgmap.array)]))

            # Write output
            deleteat!(thisgmap, .~thisgmap.array)
            CSV.write(ofile, thisgmap)
            end
        end
    end
end
end
