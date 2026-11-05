
# Load the library
using CSV
using DataFrames
using LocalAncestry

# Input
blocksizes = 0.00001:0.00005:0.002
replicates = 1:30
markersets = ["array", "sequence"]

# Constants
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"
include("src/julia/collection.jl")

# Set up the genomic map
gmap = CSV.read(gmappath, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2 = CSV.read(gmappath2, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2.array .= true
leftjoin!(gmap, gmap2, on=names(gmap))
gmap.array[ismissing.(gmap.array)] .= false

#
for r in replicates, m in markersets
    
    # Set paths
    refpath = "data/replicate_$(r)/ancestral_$(m)_10_10.vcf.gz"
    admpath = "data/replicate_$(r)/admixed_$(m)_10_10.vcf.gz"
    ancpath = "data/replicate_$(r)/global_ancestries.csv"
    trupath = "data/replicate_$(r)/local_ancestries.csv"

    # Read true local ancestries
    tru = CSV.read(trupath, DataFrame)
    tru.basepairs = UnitRange{Int}.(parse.(Int, replace.(tru.basepairs, r"\:.+" => "")), parse.(Int, replace.(tru.basepairs, r".+\:" => "")))

    @threads for b in blocksizes
        thisgmap = copy(gmap)
        est = localancestry(refpath, admpath, ancpath, threshold = Inf, maf = 0.0, maxblockfraction = b)
        accuracies = evaluate(est, tru, thisgmap)
        thisgmap.accuracy = accuracies.pcorrect
        thisgmap.incorrect = accuracies.pincorrect
        thisgmap.unassigned = accuracies.punassigned
        println(mean(thisgmap.accuracy[Bool.(thisgmap.array)]))

        # Write output
        deleteat!(thisgmap, .~thisgmap.array)
        CSV.write("results/blocks/$(m)_$(r)_$(b).csv", thisgmap)

    end
end
