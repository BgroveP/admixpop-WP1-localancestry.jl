

using CSV
using DataFrames
using LocalAncestry
using QGIO
using Base.Threads

# Get ARGS
vcfpath = ARGS[1]
truelapath = ARGS[2]
outpath = ARGS[3]

# Read custom function
include("READ_FLARE.jl")
include("DUMP.jl")

# Test
#vcfpath = "results/flare/29_array_10_10.anc.vcf.gz"
#truelapath = "data/replicate_29/local_ancestries.csv"
#outpath = "results/flare/29_array_10_10.csv"

# Hard-coded paths
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"

# Read
est = read_flare(vcfpath)
tru = CSV.read(truelapath, DataFrame)
gmap = CSV.read(gmappath, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2 = CSV.read(gmappath2, DataFrame, header = ["chromosome", "position(bp)", "distance"])

# Determine sequence and array loci
gmap2.array .= true
leftjoin!(gmap, gmap2, on=names(gmap))
gmap.array[ismissing.(gmap.array)] .= false

# Get accuracy
accuracies = evaluate(est, tru, gmap)
gmap.accuracy = accuracies.pcorrect
gmap.incorrect = accuracies.pincorrect
gmap.unassigned = accuracies.punassigned

# Write output
CSV.write(outpath, gmap)