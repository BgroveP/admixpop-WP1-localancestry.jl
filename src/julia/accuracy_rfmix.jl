

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
# VCF whose loci are used to translate window end positions (the target VCF of the run)
locipath = length(ARGS) >= 4 ? ARGS[4] : replace(truelapath, "local_ancestries.csv" => "admixed_sequence_10_10.vcf.gz")

# Read custom function
println("Reading custom functions")
include("read_rfmix.jl")
include("collection.jl")

# Test
#vcfpath = "results/localancestry/rfmix/2_sequence_6_6.msp.tsv"
#truelapath = "data/replicate_2/local_ancestries.csv"
#outpath = "results/rfmix/2_sequence_6_6.csv"

# Hard-coded paths
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"

# Read
println("Reading data")
est = read_rfmix2(vcfpath, QGIO.loci(locipath))
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
accuracies = evaluate(est, tru, gmap)
gmap.accuracy = accuracies.pcorrect
gmap.incorrect = accuracies.pincorrect
gmap.unassigned = accuracies.punassigned
println(mean(gmap.accuracy[Bool.(gmap.array)]))

# Maybe subset
arrayfile = contains(vcfpath, "_array_")
arrayfile ? deleteat!(gmap, .~gmap.array) : nothing

# Write output
CSV.write(outpath, gmap)
