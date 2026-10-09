

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
include("read_sparsepainter.jl")
include("collection.jl")

# Test
# include("src/julia/read_sparsepainter.jl")
# include("src/julia/collection.jl")
#vcfpath = "results/localancestry/sparsepainter/1_array_10_6_prob.txt.gz"
#truelapath = "data/replicate_1/local_ancestries.csv"
#outpath = "results/accuracy/sparsepainter/1_array_10_6.csv"

# Hard-coded paths
gmappath = "data/rfmix_sequence.tsv"
gmappath2 = "data/rfmix_array.tsv"

# Read
est = read_sparsepainter(vcfpath, "chr22")
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
accuracies = evaluate2(est, tru, gmap)
gmap.accuracy = accuracies.pcorrect
gmap.incorrect = accuracies.pincorrect
gmap.unassigned = accuracies.punassigned
println(mean(gmap.accuracy[Bool.(gmap.array)]))

# Write output
CSV.write(outpath, gmap)
# combine(groupby(tru, :ancestry), :basepairs => (x -> sum(sum.(x))))