# Load packages
using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFrames

# Paths 
resultspath = "results/summary/"

# Sensitivity plots
grcols = [:software, :markerset, :admgen]
transform!(groupby(olong, grcols),
	[:memory, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_memory_mean,
	[:time, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_time_mean,
	[:accuracy, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

# Sensitivity to number of target individuals
analyses.target = parse.(Int, replace.(analyses.Ntarget, "ADM10" => "2000", "ADM8-10" => "4000", "ADM6-10" => "6000"))
for stat in ["time", "memory", "accuracy"]
	tmp = combine(groupby(analyses, [:lociset, :Nreference, :target, :program]),
		"relative_$(stat)" => mean => "$(stat)")
	tmpplot = data(tmp) *
			  mapping(
				  :target => "Number of target individuals",
				  "$(stat)",
				  col = :lociset => "Marker set",
				  row = :Nreference,
				  color = :program => "LAI tool") *
			  visual(Lines)
	fig = draw(tmpplot, axis = (; xticklabelrotation = pi / 4, limits = (nothing, (minimum(tmp[:, stat]), maximum(tmp[tmp.program .!= "FLARE", stat])))))
	save("plots/article/WP1/sensitivity_$(stat)_target.pdf", fig)
end

# Sensitivity to number of reference individuals
grcols = [:program, :lociset, :Ntarget]
transform!(groupby(analyses, grcols),
	[:memory, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_memory_mean,
	[:time, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_time_mean,
	[:accuracy, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean
analyses.reference = parse.(Int, replace.(analyses.Nreference, "REF10" => "6000", "REF8-10" => "18000", "REF6-10" => "30000"))
for stat in ["time", "memory", "accuracy"]
	tmp = combine(groupby(analyses, [:lociset, :Ntarget, :reference, :program]),
		"relative_$(stat)" => mean => "$(stat)")
	tmpplot = data(tmp) *
			  mapping(
				  :reference => "Number of reference individuals",
				  "$(stat)",
				  col = :lociset => "Marker set",
				  row = :Ntarget,
				  color = :program => "LAI tool") *
			  visual(Lines)
	fig = draw(tmpplot, axis = (; xticklabelrotation = pi / 4, limits = (nothing, (minimum(tmp[:, stat]), maximum(tmp[tmp.program .!= "FLARE", stat])))))
	save("plots/article/WP1/sensitivity_$(stat)_reference.pdf", fig)
end

# Sensitivity to number of markers
grcols = [:program, :Nreference, :Ntarget]
transform!(groupby(analyses, grcols),
	[:memory, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_memory_mean,
	[:time, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_time_mean,
	[:accuracy, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean
combine(groupby(analyses, [:program, :lociset]),
	[:relative_memory, :relative_time] .=> median,
	[:relative_memory, :relative_time] .=> (x -> sqrt(var(x))))
for stat in ["time", "memory", "accuracy"]
	tmp = combine(groupby(analyses, [:lociset, :Ntarget, :reference, :program]),
		"relative_$(stat)" => mean => "$(stat)")
	tmpplot = data(tmp) *
			  mapping(
				  :lociset => "Number of reference individuals",
				  "$(stat)",
				  col = :reference => "Marker set",
				  row = :Ntarget,
				  color = :program => "LAI tool") *
			  visual(Lines)
	fig = draw(tmpplot, axis = (; xticklabelrotation = pi / 4, limits = (nothing, (minimum(tmp[:, stat]), maximum(tmp[tmp.program .!= "FLARE", stat])))))
	save("plots/article/WP1/sensitivity_$(stat)_markers.pdf", fig)
end

#### ----- Second attempt at plotting ----- #####

# Sensitivity plots
grcols = [:program, :lociset, :Nreference]
transform!(groupby(analyses, grcols),
	[:memory, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_memory_mean,
	[:time, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_time_mean,
	[:accuracy, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

# Sensitivity to number of target individuals
df = stack(analyses, [:relative_time, :relative_memory, :relative_accuracy])[:, [:program, :lociset, :Nreference, :Ntarget, :variable, :value]]
c_order = unique(sort(analyses, :time).Ntarget)
deleteat!(df, df.value .>= 5)
df.variable = replace.(df.variable, "relative_time" => "Time", "relative_memory" => "Memory", "relative_accuracy" => "Accuracy")
tmpdata = combine(groupby(df, [:Ntarget, :program, :variable]), :value => mean => :value)
tmpplot =
	data(df) *
	mapping(
		:Ntarget => sorter(c_order) => " ",
		:value => "Relative increase",
		color = :program => "LAI tool",
		col = :variable) *
	visual(RainClouds, plot_boxplots = false, show_median = false, clouds = nothing, side_nudge = 0.0) +
	(data(tmpdata) * mapping(:Ntarget => sorter(c_order), :value => "Relative increase", color = :program, col = :variable) * visual(Lines))
fig1 = draw(tmpplot, legend=(; position=:top), axis = (; xticklabelrotation = pi/4), facet=(; linkyaxes=:none))
save("plots/article/WP1/sensitivity_allstats_target.pdf", fig1)
sort(combine(groupby(df, [:program, :variable]), :value .=> [mean, minimum, maximum]), [:variable, :value_mean])

# Sensitivity to number of reference individuals
grcols = [:program, :lociset, :Ntarget]
transform!(groupby(analyses, grcols),
	[:memory, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_memory_mean,
	[:time, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_time_mean,
	[:accuracy, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

df = stack(analyses, [:relative_time, :relative_memory, :relative_accuracy])[:, [:program, :lociset, :Nreference, :Ntarget, :variable, :value]]
c_order = unique(sort(analyses, :time).Nreference)
sort(df, :value)
deleteat!(df, df.value .>= 21)
df.variable = replace.(df.variable, "relative_time" => "Time", "relative_memory" => "Memory", "relative_accuracy" => "Accuracy")
tmpdata = combine(groupby(df, [:Nreference, :program, :variable]), :value => mean => :value)
tmpplot2 =
	data(df) *
	mapping(
		:Nreference => sorter(c_order) => " ",
		:value => "Relative increase",
		color = :program => "LAI tool",
		col = :variable) *
	visual(RainClouds, plot_boxplots = false, show_median = false, clouds = nothing, side_nudge = 0.0) +
	(data(tmpdata) * mapping(:Nreference => sorter(c_order), :value => "Relative increase", color = :program, col = :variable) * visual(Lines))
fig2 = draw(tmpplot2, legend=(; position=:top), axis = (; xticklabelrotation = pi/4), facet=(; linkyaxes=:none))
save("plots/article/WP1/sensitivity_allstats_reference.pdf", fig2)
sort(combine(groupby(df, [:program, :variable]), :value .=> [mean, minimum, maximum]), [:variable, :value_mean])

# Sensitivity to number of markers
grcols = [:program, :Nreference, :Ntarget]
transform!(groupby(analyses, grcols),
	[:memory, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_memory_mean,
	[:time, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_time_mean,
	[:accuracy, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

df = stack(analyses, [:relative_time, :relative_memory, :relative_accuracy])[:, [:program, :lociset, :Nreference, :Ntarget, :variable, :value]]
c_order = unique(sort(analyses, :time).lociset)
sort(df, :value)
deleteat!(df, (df.value .>= 15) .| (df.program .== "Recomb-Mix"))
df.variable = replace.(df.variable, "relative_time" => "Time", "relative_memory" => "Memory", "relative_accuracy" => "Accuracy")
tmpdata = combine(groupby(df, [:lociset, :program, :variable]), :value => mean => :value)
tmpplot3 =
	data(df) *
	mapping(
		:lociset => sorter(c_order) => " ",
		:value => "Relative increase",
		color = :program => "LAI tool",
		col = :variable) *
	visual(RainClouds, plot_boxplots = false, show_median = false, clouds = nothing, side_nudge = 0.0) +
	(data(tmpdata) * mapping(:lociset => sorter(c_order), :value => "Relative increase", color = :program, col = :variable) * visual(Lines))
fig3 = draw(tmpplot3, legend=(; position=:top), axis = (; xticklabelrotation = pi/4), facet=(; linkyaxes=:none))
save("plots/article/WP1/sensitivity_allstats_markers.png", fig3)
sort(combine(groupby(df, :program), :value .=> [mean, minimum, maximum]), :value_mean)



#### ----- Third attempt at plotting ----- #####

# Sensitivity plots
grcols = [:program, :lociset, :Nreference]
transform!(groupby(analyses, grcols),
	[:memory, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_memory_mean,
	[:time, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_time_mean,
	[:accuracy, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(z .== "ADM10")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

# Sensitivity to number of target individuals
df = stack(analyses, [:relative_time, :relative_memory, :relative_accuracy])[:, [:program, :lociset, :Nreference, :Ntarget, :variable, :value]]
c_order = unique(sort(analyses, :time).Ntarget)
df.variable = replace.(df.variable, "relative_time" => "Time", "relative_memory" => "Memory", "relative_accuracy" => "Accuracy")
df.x = df.Ntarget
df.group .= "Target"
plotdf = copy(df)

# Sensitivity to number of reference individuals
grcols = [:program, :lociset, :Ntarget]
transform!(groupby(analyses, grcols),
	[:memory, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_memory_mean,
	[:time, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_time_mean,
	[:accuracy, :Nreference, :Ntarget] => ((x, y, z) -> mean(x[(y .== "REF10")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

df = stack(analyses, [:relative_time, :relative_memory, :relative_accuracy])[:, [:program, :lociset, :Nreference, :Ntarget, :variable, :value]]
c_order = unique(sort(analyses, :time).Nreference)
df.variable = replace.(df.variable, "relative_time" => "Time", "relative_memory" => "Memory", "relative_accuracy" => "Accuracy")

df.x = df.Nreference
df.group .= "Reference"
plotdf = append!(plotdf, df)

# Sensitivity to number of markers
grcols = [:program, :Ntarget, :Nreference]
transform!(groupby(analyses, grcols),
	[:memory, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_memory_mean,
	[:time, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_time_mean,
	[:accuracy, :lociset, :Ntarget] => ((x, y, z) -> mean(x[(y .== "Array")])) => :relative_accuracy_mean)
analyses.relative_memory = analyses.memory ./ analyses.relative_memory_mean
analyses.relative_time = analyses.time ./ analyses.relative_time_mean
analyses.relative_accuracy = analyses.accuracy ./ analyses.relative_accuracy_mean

df = stack(analyses, [:relative_time, :relative_memory, :relative_accuracy])[:, [:program, :lociset, :Nreference, :Ntarget, :variable, :value]]
c_order = unique(sort(analyses, :time).lociset)
df.variable = replace.(df.variable, "relative_time" => "Time", "relative_memory" => "Memory", "relative_accuracy" => "Accuracy")

df.x = df.lociset
df.group .= "Marker"
plotdf = append!(plotdf, df)
deleteat!(plotdf, plotdf.value .>= 15)
tmpdf = combine(groupby(plotdf, [:x, :variable, :program, :group]), :value => mean => :value)
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
sort!(plotdf, :x)
### Do the plot
tmpplot =
	data(plotdf) *
	mapping(
		:x => sorter(c_order) => "",
		:value => "Relative increase",
		color = :program => "LAI tool",
		row = :variable,
		col = :group) *
	visual(RainClouds, plot_boxplots = false, show_median = false, clouds = nothing, side_nudge = 0.0) +
	(data(tmpdf) * mapping(
		:x => sorter(c_order) => "",
		:value => "Relative increase",
		color = :program => "LAI tool",
		row = :variable,
		col = :group) * visual(Lines))
fig = draw(tmpplot, 
legend=(; position=:top), 
axis = (; xticklabelrotation = pi/4), 
facet=(; linkxaxes=:colwise, linkyaxes=:rowwise))


xlims!(fig.figure.content[1], 6.8,8.2)
xlims!(fig.figure.content[2], 6.8,8.2)
xlims!(fig.figure.content[3], 6.8,8.2)
xlims!(fig.figure.content[4], 3.8,6.2)
xlims!(fig.figure.content[5], 3.8,6.2)
xlims!(fig.figure.content[6], 3.8,6.2)
xlims!(fig.figure.content[7], 0.8,3.2)
xlims!(fig.figure.content[8], 0.8,3.2)
xlims!(fig.figure.content[9], 0.8,3.2)
save("plots/article/WP1/sensitivity_allstats_allgroups2.pdf", fig)

# Scatter plot
tmp = analyses[:,[:replicate, :referenceset, :program, :lociset, :Nreference, :Ntarget, :accuracy]]
tmp = stack(unstack(tmp, :program, :accuracy), [:RFMix, :FLARE, Symbol("Recomb-Mix"), Symbol("LocalAncestry.jl")])
tmp = sort(vcat(tmp, tmp, tmp, tmp), [:replicate, :referenceset, :lociset, :Nreference, :Ntarget, :variable])

tmp2 = copy(tmp)
tmp2.value2 = tmp2.value
tmp2.variable2 = tmp2.variable
select!(tmp2, [:replicate, :referenceset, :lociset, :Nreference, :Ntarget, :variable2, :value2])

tmp.variable2 = repeat(["FLARE", "LocalAncestry.jl", "RFMix", "Recomb-Mix"], outer = Int(nrow(tmp)/4))
tmp = sort(innerjoin(tmp, tmp2, on = [:replicate, :referenceset, :lociset, :Nreference, :Ntarget, :variable2]))
deleteat!(tmp, ismissing.(tmp.value .* tmp.value2))
tmp.programs = "1: " .* tmp.variable .* ". 2: " .* tmp.variable2
tmp.program = tmp.variable
tmpplot =
	data(tmp) *
	mapping(
		:value => "Accuracy",
		:value2 => "Accuracy",
		color = :program => "LAI tool",
		layout  = :variable2) *
	visual(Scatter) +
	data(DataFrame(value = [minimum(tmp.value), maximum(tmp.value)], value2 = [minimum(tmp.value), maximum(tmp.value)])) * 
	mapping(
		:value => "Accuracy",
		:value2 => "Accuracy") * visual(Lines)

fig = draw(tmpplot, legend=(; position=:top), 
axis = (; xticklabelrotation = pi/4), 
facet=(; linkxaxes=:colwise, linkyaxes=:rowwise))
save("plots/article/WP1/accuracy_allstats_allgroups.pdf", fig)

### Create only the subplot with RFmix as x-axis
tmp = analyses[:,[:replicate, :referenceset, :program, :lociset, :Nreference, :Ntarget, :accuracy, :accuracy_hapmap3]]
tmp.accuracy = [l == "Sequence" ? tmp.accuracy[i] : tmp.accuracy_hapmap3[i]  for (i, l) in enumerate(tmp.lociset)] 
select!(tmp, Not(:accuracy_hapmap3))
tmp = stack(unstack(tmp, :program, :accuracy), [:RFMix, :FLARE, Symbol("Recomb-Mix"), Symbol("LocalAncestry.jl")])
tmp = sort(tmp, [:replicate, :referenceset, :lociset, :Nreference, :Ntarget, :variable])

tmp2 = copy(tmp)
tmp2.value2 = tmp2.value
tmp2.variable2 = tmp2.variable
select!(tmp2, [:replicate, :referenceset, :lociset, :Nreference, :Ntarget, :variable2, :value2])

tmp.variable2 = repeat(["RFMix"], outer = nrow(tmp))
tmp = sort(innerjoin(tmp, tmp2, on = [:replicate, :referenceset, :lociset, :Nreference, :Ntarget, :variable2]))
deleteat!(tmp, ismissing.(tmp.value .* tmp.value2))
tmp.program = tmp.variable
tmp2 = copy(tmp)
tmp2 = deleteat!(tmp2, (tmp2.variable2 .!= "RFMix") )
tmp2 = deleteat!(tmp2, (tmp2.variable .== "RFMix"))
tmpplot =
	data(tmp2) *
	mapping(
		:value2 => "Accuracy (RFMix)",
		:value => "Accuracy (Color)",
		color = :variable => "LAI tool",
		marker = :lociset => "Marker set") *
	visual(Scatter) +
	data(DataFrame(value = [minimum(tmp2.value2), maximum(tmp2.value2)], value2 = [minimum(tmp2.value2), maximum(tmp2.value2)])) * 
	mapping(
		:value2 => "Accuracy (RFMix)",
		:value => "Accuracy (Color)") * visual(Lines)

fig = draw(tmpplot, legend=(; position=:top), 
axis = (; xticklabelrotation = pi/4), 
facet=(; linkxaxes=:colwise, linkyaxes=:rowwise))

save("plots/article/WP1/accuracy_allstats_allgroups2.pdf", fig)

# Numbers for the article
### Average accuracies across all scenarios
combine(groupby(analyses, :program), 
[:accuracy, :accuracy_hapmap3, :lociset] => ((x,y, z) -> mean([x[z .== "Sequence"]; y])) => "mean",
[:accuracy, :accuracy_hapmap3, :lociset] => ((x,y, z) -> std([x[z .== "Sequence"]; y])) => "std",
) |> sort

### Average accuracy for array -> array
combine(groupby(analyses[analyses.lociset .== "Array", :], :program), :accuracy_hapmap3 .=> [mean, std]) |> sort

### Average accuracy for sequence -> sequence
combine(groupby(analyses[analyses.lociset .== "Sequence", :], :program), :accuracy .=> [mean, std]) |> sort

### Average accuracy for sequence -> sequence
combine(groupby(analyses[analyses.lociset .== "Sequence", :], :program), :accuracy_hapmap3 .=> [mean, std]) |> sort

### Average accuracy for sequence -> sequence
tmp = unstack(analyses[:, [1:5; 7:7; 12]], :lociset, :accuracy_hapmap3)
tmp.dif = tmp.Sequence -tmp.Array 
deleteat!(tmp, ismissing.(tmp.dif))
combine(groupby(tmp, :program), :dif .=> [mean, std]) |> sort

### Average time use for all scenarios - (ADM10, ADM8-10) x array->array x REF10
tmp = copy(analyses)
deleteat!(tmp, ((tmp.lociset .== "Array") .& (tmp.Nreference .== "REF10")) .& ((tmp.Ntarget .== "ADM10") .| (tmp.Ntarget .== "ADM8-10")))
combine(groupby(tmp, [:program]), :time .=> [mean, std])


### Average time use for all scenarios - (ADM10, ADM8-10) x array->array x REF10
tmp = copy(analyses)
tmp2 = combine(groupby(tmp, [:program, :lociset, :Ntarget, :Nreference]), :memory .=> [mean, std])
sort(tmp2[tmp2.Ntarget .== "ADM10",:], :memory_mean)
sort(tmp2[(tmp2.Ntarget .!= "ADM10") .& (tmp2.lociset .!= "Array"),:], :memory_mean)


### Time comparison
tmp = copy(analyses)
tmp2 = combine(groupby(tmp[tmp.lociset .== "Sequence",:], [:program]), :time .=> [mean, std])

### Memory comparison
tmp = copy(analyses)
tmp2 = combine(groupby(tmp[tmp.lociset .== "Sequence",:], [:program]), :memory .=> [mean, std])
