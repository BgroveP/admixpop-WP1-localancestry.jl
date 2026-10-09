

using CSV
using DataFrames
using LocalAncestry
using QGIO
using Base.Threads
using StatsBase

# test
# ARGS = ["results/localancestry/localancestryjl/1_sequence_10_10_0.0.csv",
#        "data/replicate_1/local_ancestries.csv",
#        "results/accuracy/localancestryjl/1_sequence_10_10_0.0.csv"]

# Get ARGS
vcfpath = ARGS[1]
truelapath = ARGS[2]
outpath = ARGS[3]


# Hard-coded variables
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"
include("collection.jl")

# Read and calculate accuracy
est = CSV.read(vcfpath, DataFrame)
tru = CSV.read(truelapath, DataFrame)
gmap = CSV.read(gmappath, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2 = CSV.read(gmappath2, DataFrame, header = ["chromosome", "position(bp)", "distance"])

# Convert basepairs strings to ranges
est.basepairs = UnitRange{Int}.(parse.(Int, replace.(est.basepairs, r"\:.+" => "")), parse.(Int, replace.(est.basepairs, r".+\:" => "")))
tru.basepairs = UnitRange{Int}.(parse.(Int, replace.(tru.basepairs, r"\:.+" => "")), parse.(Int, replace.(tru.basepairs, r".+\:" => "")))

# Determine sequence and array loci
gmap2.array .= true
leftjoin!(gmap, gmap2, on=names(gmap))
gmap.array[ismissing.(gmap.array)] .= false

# Get accuracy
accuracies = evaluate2(est, tru, gmap)
gmap.accuracy = accuracies.pcorrect
gmap.incorrect = accuracies.pincorrect
gmap.unassigned = accuracies.punassigned
println(mean(gmap.accuracy[Bool.(gmap.array)]))

# Write output
CSV.write(outpath, gmap)