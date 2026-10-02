

using CSV
using DataFrames
using LocalAncestry
using QGIO
using Base.Threads
using StatsBase

# Get ARGS
vcfpath = ARGS[1]
truelapath = ARGS[2]
outpath = ARGS[3]

# Read custom function
include("read_flare.jl")
include("collection.jl")

# Test
#vcfpath = "results/localancestry/flare/1_array_6_6.anc.vcf.gz"
#truelapath = "data/replicate_1/local_ancestries.csv"
#outpath = "results/accuracy/flare/1_array_6_6.csv"

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

# Convert basepairs strings to ranges
tru.basepairs = UnitRange{Int}.(parse.(Int, replace.(tru.basepairs, r"\:.+" => "")), parse.(Int, replace.(tru.basepairs, r".+\:" => "")))

# Get accuracy
accuracies = evaluate(est, tru, gmap)
gmap.accuracy = accuracies.pcorrect
gmap.incorrect = accuracies.pincorrect
gmap.unassigned = accuracies.punassigned
println(mean(gmap.accuracy[Bool.(gmap.array)]))

# Write output
CSV.write(outpath, gmap)

