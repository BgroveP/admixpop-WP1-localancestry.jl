
# Loac packages
using CSV
using DataFrames
using StatsBase
using Plots

# Input
defaultt = 0.6

# 
softwares = readdir("results/accuracy")
softwares = softwares[softwares .!= "recombmix"]
nsoftwares = length(softwares)

# Read
o = Vector{DataFrame}(undef, nsoftwares)
for (i,s) in enumerate(softwares)
    println(s)
    b = readdir("results/accuracy/$(s)", join = true)
    o[i] = CSV.read(b, DataFrame)
    o[i].file = repeat(replace.(b, "results/accuracy/" =>""), inner = Int(nrow(o[i])/length(b)))
end

# Collect
o = vcat(o...)

# Get information from file name
parts = split.(basename.(o.file), [['_', '.']])
o.software  = dirname.(o.file)
o.replicate = parse.(Int, getindex.(parts, 1))
o.markerset    = getindex.(parts, 2)
o.ancgen    = parse.(Int, getindex.(parts, 3))
o.admgen    = parse.(Int, getindex.(parts, 4))
o.threshold  .= defaultt
o.threshold[o.software .== "localancestryjl"]    = parse.(Int, getindex.(parts[o.software .== "localancestryjl"], 5)) + parse.(Int, getindex.(parts[o.software .== "localancestryjl"], 6))/10 

# Get overall accuracies
o.keep = (o.array .& (o.markerset .== "array")) .| (o.markerset .== "sequence")
deleteat!(o, .~o.keep)
o2 = combine(groupby(o, ["software", "ancgen", "admgen", "markerset", "threshold"]), :accuracy .=> [mean, minimum, maximum])
CSV.write("results/accuracy/summary.csv", o2)
#unstack(o2[:,Not([:accuracy_minimum,:accuracy_maximum])], :software, :accuracy_mean)

# Plot accuracy across genome (array)
tplot = scatter(o2[(o2.software .== "flare") .& (o2.markerset .== "array"), "position(bp)"], o2.accuracy_mean[(o2.software .== "flare")  .& (o2.markerset .== "array")], label = "flare")
scatter!(tplot, o2[(o2.software .== "localancestryjl") .& (o2.markerset .== "array"), "position(bp)"], o2.accuracy_mean[(o2.software .== "localancestryjl")  .& (o2.markerset .== "array")], label = "localancestryjl")
scatter!(tplot, o2[(o2.software .== "rfmix") .& (o2.markerset .== "array"), "position(bp)"], o2.accuracy_mean[(o2.software .== "rfmix")  .& (o2.markerset .== "array")], label = "rfmix")
#plot!(tplot, o2[(o2.software .== "recombmix") .& (o2.markerset .== "array"), "position(bp)"], o2.accuracy_mean[(o2.software .== "recombmix")  .& (o2.markerset .== "array")], label = "recombmix")
#plot!(tplot, o2[(o2.software .== "sparsepainter") .& (o2.markerset .== "array"), "position(bp)"], o2.accuracy_mean[(o2.software .== "sparsepainter")  .& (o2.markerset .== "array")], label = "sparsepainter")
savefig(tplot, "plots/Accuracy_chromosomewide.png")

# Plot across threshold for localancestryjl
olajl = deleteat!(o, o.software .!= "localancestryjl")
olajl.position = olajl[:, "position(bp)"]
olajl2 = unstack(select(olajl, [:position, :markerset, :accuracy, :ancgen, :admgen, :replicate, :threshold]), :threshold, :accuracy)
olajl3 = combine(groupby(olajl2, [:position, :markerset, :ancgen, :admgen]), ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"] .=> mean .=> "t" .* ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"])

cols = names(olajl3, 5:ncol(olajl3))
tplot = scatter(olajl3.position, Matrix(olajl3[:, 5:end]), label = permutedims(cols))
savefig(tplot, "plots/Accuracy_chromosomewide_localancestryjl_thresholds.png");

