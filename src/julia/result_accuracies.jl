
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
softwares = ["flare", "recombmix", "gnomix","localancestryjl","rfmix","sparsepainter"]
softwares = softwares[isdir.("results/accuracy/" .* softwares)]
softwares = softwares[softwares .!= "flareclust"]
nsoftwares = length(softwares)

# Read accuracies
o = Vector{Vector{DataFrame}}(undef, nsoftwares)
@threads for i in axes(softwares, 1)
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
olong_sds = combine(groupby(olong_means2, ["software", "ancgen", "admgen", "markerset", "array", "threshold"]), :accuracy => (x -> sqrt(var(x))) => :sd)
owide_sds = sort(unstack(olong_sds, :software, :sd), [:threshold, :array, :ancgen, :admgen, :markerset])
CSV.write("results/summary/accuracy_sds.csv", owide_sds)

mean(olong_sds.sd)*100
minimum(olong_sds.sd)*100
maximum(olong_sds.sd)*100

# Write table
tbl_ancgen = 10
meansandsds = leftjoin(olong_means, olong_sds, on = names(olong_means)[1:6])
deleteat!(meansandsds, meansandsds.ancgen .!= tbl_ancgen)
meansandsds.admgen = string.(meansandsds.admgen) 

meansandsdswider = unstack(meansandsds[:, [:software, :markerset, :accuracy, :array, :admgen, :threshold]], :admgen, :accuracy)
sort!(meansandsdswider, [:markerset, :software])
colorder = ["10", "8", "6"]
meansandsdswider = meansandsdswider[:, [["software", "markerset", "threshold", "array"]; colorder]]
meansandsdswider.outstring .= ""
for i in axes(meansandsdswider, 1)
    t = [meansandsdswider[i, j] for j in 5:(ncol(meansandsdswider) - 1)]
    meansandsdswider[i,end] = join(round.(100*t, digits = 1), " & ")
end
sort!(meansandsdswider, [:markerset, :array, :software, :threshold])
CSV.write("tbloutstring_accuracy.csv", meansandsdswider[:, [:threshold, :array, :software, :markerset, :outstring]])

# Plot histogram of accuracies
olong.software = replace.(string.(olong.software), 
							"flare" => "FLARE", "localancestryjl" => "LocalAncestry.jl", 
							"gnomix" => "Gnomix", "recombmix" => "Recomb-Mix", 
							"rfmix" => "RFMix", "sparsepainter" => "SparsePainter")
olong.markerset = replace.(string.(olong.markerset), "array" => "Array", "sequence" => "Sequence")

tmpplot =
	data(olong) *
	mapping(
		:accuracy => "Accuracy",
        group = :markerset,
        row = :software,
		color = :software => "LAI tool") *
	visual(Density) 
fig = draw(tmpplot)

save("plots/accuracy_density.png", fig)
save("plots/accuracy_density.pdf", fig)

# Plot
obp = combine(groupby(olong, ["position(bp)", "software", "array", "markerset", "ancgen", "admgen"]), :accuracy => mean => :accuracy)
obp.x = obp[!,"position(bp)"]
sort!(obp, :x)

deleteat!(obp, obp.ancgen .!= 10)
deleteat!(obp, obp.admgen .!= 10)
deleteat!(obp, .~obp.array)

tmpplot =
	data(obp) *
	mapping(
        :x,
		:accuracy => "Accuracy",
        col = :markerset,
        color = :software => "LAI tool") *
	visual(Lines) 
fig = draw(tmpplot)

save("plots/accuracy_acrossgenome.png", fig)
save("plots/accuracy_acrossgenome.pdf", fig)
