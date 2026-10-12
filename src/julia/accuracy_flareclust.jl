

using CSV
using DataFrames
using LocalAncestry
using QGIO
using Base.Threads
using StatsBase
using Combinatorics

# Get ARGS
vcfpath = ARGS[1]
truelapath = ARGS[2]
outpath = ARGS[3]

# Read custom function
include("read_flare.jl")
include("collection.jl")

# Test
#include("src/julia/read_flare.jl")
#include("src/julia/collection.jl")
#vcfpath = "results/localancestry/flareclust/1_array_6_6.anc.vcf.gz"
#truelapath = "data/replicate_1/local_ancestries.csv"
#outpath = "results/accuracy/flareclust/1_array_6_6.csv"

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

# Translate ancestries
est.ancestry = string.(est.ancestry)
M, lab1, lab2 = overlap_matrix(tru, est)

best = argmax(permutations(1:3)) do p
    sum(M[i, p[i]] for i in 1:3)
end

popdict = Dict(lab2 .=> lab1[best])
for i in axes(est, 1)
    est.ancestry[i] = popdict[est.ancestry[i]]
end

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
