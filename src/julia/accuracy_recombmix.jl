

using CSV
using DataFrames
using LocalAncestry
using QGIO
using Base.Threads

include("READ_RECOMBMIX.jl")


# Get ARGS
vcfpath = ARGS[1]
truelapath = ARGS[2]
outpath = ARGS[3]

# Hard-coded variables
chromosome = "chr22"
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"

# Transform to supported format
read_recombmix(vcfpath, outpath, chromosome)

# Read and calculate accuracy
est = CSV.read(outpath, DataFrame)
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