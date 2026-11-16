
# Load packages
using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFrames
using StatsBase

# Usage: julia --project=. src/julia/result_sensitivity.jl <tool1,tool2,...> <default threshold>

# Input
resultspath = "results/summary/"
softwares = length(ARGS) >= 1 ? String.(split(ARGS[1], ",")) :
            ["flare", "recombmix", "gnomix", "localancestryjl", "rfmix", "sparsepainter"]
defaultt = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.66
mkpath("plots")
tooldict = Dict("flare" => "FLARE", "flareclust" => "FLARE (clustered)", "localancestryjl" => "LocalAncestry.jl",
                "gnomix" => "Gnomix", "recombmix" => "Recomb-Mix",
                "rfmix" => "RFMix", "sparsepainter" => "SparsePainter")

# Read
accuracy = CSV.read(resultspath * "accuracy.csv", DataFrame)
RAMuse = CSV.read(resultspath * "RAM_use_mean.csv", DataFrame)
timeuse = CSV.read(resultspath * "time_use_mean.csv", DataFrame)

# Add statistic column
accuracy.statistic .= "accuracy"
RAMuse.statistic .= "RAM"
timeuse.statistic .= "time"

# Remove the Sequence -> Array rows
deleteat!(accuracy, accuracy.array .& (accuracy.markerset .== "sequence"))
softwares = intersect(softwares, names(accuracy), names(RAMuse), names(timeuse))
keepcols = ["ancgen", "admgen", "markerset", "threshold", "statistic"; softwares]
o = vcat(accuracy[:, keepcols], RAMuse[:, keepcols], timeuse[:, keepcols])

# Remove non-default thresholds
deleteat!(o, o.threshold .!= defaultt)

# Stack
olong = dropmissing(stack(o, softwares, variable_name = "software"), :value)

# Calculate the sensitivity to reference panel size
grcols = [:software, :markerset, :admgen, :statistic]
transform!(groupby(olong, grcols), [:value, :ancgen] => ((x, y) -> mean(x[(y .== 10)])) => :relanc)
olong.relanc = olong.value ./ olong.relanc

# Calculate the sensitivity to the number of markers
grcols = [:software, :ancgen, :admgen, :statistic]
transform!(groupby(olong, grcols), [:value, :markerset] => ((x, y) -> mean(x[(y .== "array")])) => :relmarker)
olong.relmarker = olong.value ./ olong.relmarker

# Calculate the sensitivity to target panel size
grcols = [:software, :ancgen, :markerset, :statistic]
transform!(groupby(olong, grcols), [:value, :admgen] => ((x, y) -> mean(x[(y .== 10)])) => :reladm)
olong.reladm = olong.value ./ olong.reladm

# Rename strings
oloong = stack(olong[:, Not("value")], [:relanc, :relmarker, :reladm])
oloong.statistic = replace.(oloong.statistic, "time" => "Time", "RAM" => "Memory", "accuracy" => "Accuracy")
oloong.variable = replace.(oloong.variable, "relanc" => "Reference", "reladm" => "Target", "relmarker" => "Marker")
oloong.ancgen = replace.(string.(oloong.ancgen), "6" => "REF6-10", "8" => "REF8-10", "10" => "REF10")
oloong.admgen = replace.(string.(oloong.admgen), "6" => "ADM6-10", "8" => "ADM8-10", "10" => "ADM10")
oloong.markerset = replace.(string.(oloong.markerset), "array" => "Array", "sequence" => "Sequence")
oloong.software = [get(tooldict, s, s) for s in string.(oloong.software)]

# Make the x columns
oloong.x .= ""
for i in axes(oloong, 1)
	if oloong.variable[i] == "Reference"
		oloong.x[i] = oloong.ancgen[i]
	elseif oloong.variable[i] == "Target"
		oloong.x[i] = oloong.admgen[i]
	else
		oloong.x[i] = oloong.markerset[i]
	end
end

# Make the scatterplot
tmpdf = combine(groupby(oloong, [:variable, :software, :statistic, :x]), :value => mean => :value)
c_order = [
	 "ADM10",
	 "ADM8-10",
	 "ADM6-10",
	 "REF10",
	 "REF8-10",
	 "REF6-10",
	 "Array",
 "Sequence"
]
sort!(oloong, :x)

### Do the plot
set_theme!(theme_minimal(); fontsize = 16)   
tmpplot =
	data(oloong) *
	mapping(
		:x => sorter(c_order) => "",
		:value => "Relative increase",
		color = :software => "LAI tool",
		row = :statistic,
		col = :variable) *
	visual(RainClouds, plot_boxplots = false, show_median = false, clouds = nothing, side_nudge = 0.0) +
	(data(tmpdf) * mapping(
		:x => sorter(c_order) => "",
		:value => "Relative increase",
		color = :software => "LAI tool",
		row = :statistic,
		col = :variable) * visual(Lines))
fig = draw(tmpplot, 
    figure = (; size = (1000, 450)),
	legend=(; position = :right, tellwidth = true), axis = (; xticklabelrotation = pi/4, yscale = log10), facet=(; linkxaxes=:colwise, linkyaxes=:rowwise))


xlims!(fig.figure.content[1], 6.5,8.5)
xlims!(fig.figure.content[2], 6.5,8.5)
xlims!(fig.figure.content[3], 6.5,8.5)
xlims!(fig.figure.content[4], 3.8,6.2)
xlims!(fig.figure.content[5], 3.8,6.2)
xlims!(fig.figure.content[6], 3.8,6.2)
xlims!(fig.figure.content[7], 0.8,3.2)
xlims!(fig.figure.content[8], 0.8,3.2)
xlims!(fig.figure.content[9], 0.8,3.2)
save("plots/sensitivities.png", fig, px_per_unit = 2)
save("plots/sensitivities.pdf", fig, px_per_unit = 2)

### Accuracy of each tool against RFMix (same scenario), with the identity line
if "rfmix" in softwares && length(softwares) > 1
	others = setdiff(softwares, ["rfmix"])
	orfmix = stack(o, others, variable_name = "software", value_name = "accuracy")
	deleteat!(orfmix, orfmix.statistic .!= "accuracy")
	dropmissing!(orfmix, [:accuracy, :rfmix])
	orfmix.software = [get(tooldict, s, s) for s in string.(orfmix.software)]
	lims = extrema([orfmix.rfmix; orfmix.accuracy])

	tmpplot =
		data(orfmix) *
		mapping(
			:rfmix => "Accuracy (RFMix)",
			:accuracy => "Accuracy (Other)",
			color = :software => "LAI tool") *
		visual(Scatter) +
		data(DataFrame(x = collect(lims), y = collect(lims))) *
		mapping(:x => "Accuracy (RFMix)", :y => "Accuracy (Other)") * visual(Lines)

	fig = draw(tmpplot)
	save("plots/accuracy_vs_rfmix.png", fig)
	save("plots/accuracy_vs_rfmix.pdf", fig)
end

# Numbers for the article
### Average accuracies across all scenarios
tmp = combine(groupby(olong, [:statistic, :markerset, :software]), 
:value => mean => "mean",
:value => (x -> sqrt(var(x))) => "std"
) |> sort
tmp.mean = round.(tmp.mean, digits = 2)
tmp.std = round.(tmp.std, digits = 2)
tmp

### Average increase from using Sequence map
deleteat!(tmp, tmp.statistic .!= "accuracy")
tmpw = unstack(tmp[:, Not(:std)], :markerset, :mean)
all(in(names(tmpw)), ["sequence", "array"]) && println("Mean accuracy increase sequence vs array: ", mean(tmpw.sequence - tmpw.array))

### Increase in accuracy from using larger reference panel
mean(olong[(olong.statistic .== "accuracy") .& (olong.ancgen .!= 10),:value]) -
mean(olong[(olong.statistic .== "accuracy") .& (olong.ancgen .== 10),:value])

## Time uses with different target panels
combine(groupby(olong, [:admgen, :ancgen, :statistic]), :value .=> [mean, minimum, maximum])

## Time sensitivities to reference panels
tmpacc = combine(groupby(oloong, [:software, :statistic, :variable]), :value .=> [mean, minimum, maximum]) |> (x -> sort(x, :value_mean))
deleteat!(tmpacc, tmpacc.statistic .!= "Time")
deleteat!(tmpacc, tmpacc.variable .!= "Reference")

## Time sensitivities to target panels
tmpacc = combine(groupby(oloong, [:software, :statistic, :variable]), :value .=> [mean, minimum, maximum]) |> (x -> sort(x, :value_mean))
deleteat!(tmpacc, tmpacc.statistic .!= "Time")
deleteat!(tmpacc, tmpacc.variable .!= "Target")

# Memory uses across everything
combine(groupby(olong, :statistic), :value .=> [mean, minimum, maximum])

# Memory uses across software
tmpmem = combine(groupby(olong, [:statistic, :software]), :value .=> [mean, minimum, maximum]) |> sort
deleteat!(tmpmem, tmpmem.statistic .!= "RAM")
sort!(tmpmem, :value_mean)
