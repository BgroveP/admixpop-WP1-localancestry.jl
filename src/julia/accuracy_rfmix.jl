

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
println("
    read commandline arguments:
        - path: $(vcfpath)
        - true: $(truelapath)
        - out:  $(outpath)

    Also:
        - threads: $(nthreads())
")

# Read custom function
println("Reading custom functions")
include("read_rfmix.jl")
include("collection.jl")

# Test
#vcfpath = "results/flare/29_array_10_10.anc.vcf.gz"
#truelapath = "data/replicate_29/local_ancestries.csv"
#outpath = "results/flare/29_array_10_10.csv"

# Hard-coded paths
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"

# Read
println("Reading data")
est = read_rfmix(vcfpath, QGIO.loci(replace(truelapath, "local_ancestries.csv" => "admixed_sequence_10_10.vcf.gz")))
tru = CSV.read(truelapath, DataFrame)
gmap = CSV.read(gmappath, DataFrame, header = ["chromosome", "position(bp)", "distance"])
gmap2 = CSV.read(gmappath2, DataFrame, header = ["chromosome", "position(bp)", "distance"])

# Determine sequence and array loci
gmap2.array .= true
leftjoin!(gmap, gmap2, on=names(gmap))
gmap.array[ismissing.(gmap.array)] .= false

# Create range from string
tru.basepairs = UnitRange{Int}.(parse.(Int, replace.(tru.basepairs, r"\:.+" => "")), parse.(Int, replace.(tru.basepairs, r".+\:" => "")))

# Get accuracy
println("Calculating the accuracies")
accuracies = evaluate(est, tru, gmap)
gmap.accuracy = accuracies.pcorrect
gmap.incorrect = accuracies.pincorrect
gmap.unassigned = accuracies.punassigned
println(mean(gmap.accuracy[Bool.(gmap.array)]))

# Write output
println("Writing the output")
CSV.write(outpath, gmap)

