
# Load packages
using AlgebraOfGraphics
using Base.Threads
using CairoMakie
using CSV
using DataFrames
using Statistics
using StatsBase

# Input
defaultt = 0.6

# Get list of software
softwares = ["flare","recombmix", "gnomix","localancestryjl","rfmix","sparsepainter"]
softwares = softwares[isdir.("results/accuracy/" .* softwares)]
softwares = softwares[softwares .!= "flareclust"]
nsoftwares = length(softwares)

# Read accuracies
o = Vector{Vector{DataFrame}}(undef, nsoftwares)
for i in axes(softwares, 1)
    s = softwares[i]
    b = readdir("results/accuracy/$(s)", join = true)
    o[i] = Vector{DataFrame}(undef, length(b))
    
    for (ji, j) in enumerate(b)
        println("$(s), $(j)")
            o[i][ji] = CSV.read(j, 
            DataFrame, 
            types = [String, Int, Float64, Bool, Float64, Float64, Float64],
            stripwhitespace = true,
            strict = true)
            n = nrow(o[i][ji])
            parts = split(j, ['_', '.', '/'])

            # Get information
            o[i][ji].software  .= s
            o[i][ji].replicate .= parse(Int, parts[4])
            o[i][ji].markerset .= parts[5]
            o[i][ji].ancgen    .= parse(Int, parts[6])
            o[i][ji].admgen    .= parse(Int, parts[7])
            o[i][ji].threshold .= defaultt

            if s == "localancestryjl"
            o[i][ji].threshold .= parse(Int, parts[8]) + parse(Int, parts[9]) /10
            end
    end
end

# Collect into one dataframe
olong = reduce(vcat, vcat(o...))

# Omit locus that is never estimated (50799821) and the one never estimated for rfmix and gnomix (50780578)
olong.keep = (olong[!,"position(bp)"] .!= 50799821) .& (olong[!,"position(bp)"] .!= 50780578)
olong.keep = olong.keep .& ((olong.array .& (olong.markerset .== "array")) .| (olong.markerset .== "sequence"))
deleteat!(olong, .~olong.keep)

# Investigate whether all analyses are finished
oamount = combine(groupby(olong, ["software", "ancgen", "admgen", "markerset", "threshold", "replicate"]), :accuracy => (x -> true) => :n) 
oamountseen = stack(unstack(oamount, :software, :n), softwares, variable_name = :software, value_name = :seen)
deleteat!(oamountseen, .~ismissing.(oamountseen.seen))
deleteat!(oamountseen, (oamountseen.threshold .!= defaultt) .& (oamountseen.software .!= "localancestryjl"))
CSV.write("results/summary/missing_accuracy.csv", oamountseen)

# Get overall accuracies
olong_means = combine(groupby(olong, ["software", "ancgen", "admgen", "markerset", "array", "threshold"]), :accuracy => mean => :accuracy)
owide_means = sort(unstack(olong_means, :software, :accuracy), [:threshold, :array, :ancgen, :admgen, :markerset])
CSV.write("results/summary/accuracy.csv", owide_means)

olong_means2 = combine(groupby(olong, ["software", "ancgen", "admgen", "markerset", "array", "threshold", "replicate"]), :accuracy => mean => :accuracy)
olong_sds = combine(groupby(olong_means2, ["software", "ancgen", "admgen", "markerset", "array", "threshold"]), :accuracy => (x -> sqrt(var(x))) => :accuracy)
owide_sds = sort(unstack(olong_sds, :software, :accuracy), [:threshold, :array, :ancgen, :admgen, :markerset])
CSV.write("results/summary/accuracy_sds.csv", owide_sds)

# Plot accuracy across genome (sequence)
omostinfo = olong[(olong.ancgen .== 6) .&(olong.admgen .== 6) .& (olong.markerset .== "sequence") .& (olong.threshold .== defaultt),:]
locus_means = combine(groupby(omostinfo, ["software", "position(bp)"]), :accuracy => mean => :accuracy)
rename!(locus_means, [:software, :position, :accuracy])


tplot = data(locus_means) *
		  mapping(
			  :position => "Chr22:basepair",
              :accuracy => "Accuracy",
			  color = :software => "Software") *
		  visual(Lines) |> (x -> draw(x, axis = (; xticklabelrotation = pi / 4)))
save("plots/Accuracy_chromosomewide.png", tplot)

# Plot across threshold for localancestryjl
othreshold = olong[(olong.ancgen .== 6) .&(olong.admgen .== 6) .& (olong.markerset .== "sequence") .& (.~(olong.array)) .& (olong.software .== "localancestryjl"),:]
locus_means_threshold = combine(groupby(othreshold, ["threshold", "position(bp)"]), :accuracy => mean => :accuracy)
rename!(locus_means_threshold, [:threshold, :position, :accuracy])

tplot = data(locus_means_threshold) *
		  mapping(
			  :position => "Chr22:basepair",
              :accuracy => "Accuracy",
              group = :threshold,
			  color = :threshold => "Threshold") *
		  visual(Lines) |> (x -> draw(x, axis = (; xticklabelrotation = pi / 4)))
save("plots/Accuracy_chromosomewide_localancestryjl_thresholds.png", tplot)
