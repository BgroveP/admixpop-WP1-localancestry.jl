
# Accuracy plots: effect of reference panel size, and comparison across tools.
# Usage: julia --project=. src/julia/result_plot.jl <tool1,tool2,...> <default threshold>

using CSV
using DataFrames
using Statistics
using StatsBase
using AlgebraOfGraphics
using CairoMakie

# Input
softwares = length(ARGS) >= 1 ? String.(split(ARGS[1], ",")) :
            ["flare", "recombmix", "gnomix", "localancestryjl", "rfmix", "sparsepainter"]
default_threshold = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.66
default_ancgen = 10
mkpath("plots")

# Translations
msadict = Dict(
    ["sequence", false] => "Sequence -> Sequence",
    ["sequence", true] => "Sequence -> Array",
    ["array", true] => "Array -> Array",
)

# Read the accuracies (only tools present in the summary)
acc = CSV.read("results/summary/accuracy.csv", DataFrame, stripwhitespace = true, missingstring = ["NaN", ""])
accsd = CSV.read("results/summary/accuracy_sds.csv", DataFrame, stripwhitespace = true, missingstring = ["NaN", ""])
softwares = intersect(softwares, names(acc))
accuracies = dropmissing(stack(acc, softwares, variable_name = :software, value_name = :mean), :mean)
accuraciessd = dropmissing(stack(accsd, intersect(softwares, names(accsd)), variable_name = :software, value_name = :sd), :sd)
leftjoin!(accuracies, accuraciessd, on = [:threshold, :array, :ancgen, :admgen, :markerset, :software])

# Translate columns
accuracies.markerset = String.(accuracies.markerset)
for i in axes(accuracies, 1)
    accuracies.markerset[i] = get(msadict, [accuracies.markerset[i], accuracies.array[i]], accuracies.markerset[i])
end
accuracies.group = string.(accuracies.software, "|", accuracies.markerset, "|", accuracies.admgen)

# Accuracy vs. reference panel size
accuracy_default = accuracies[accuracies.threshold .== default_threshold, :]
nrow(accuracy_default) == 0 && error("No accuracies found for threshold $(default_threshold)")
plt = data(accuracy_default) *
      mapping(:ancgen,
              :mean => "Avg. accuracy",
              col = :admgen,
              color = :software,
              group = :group) *
      visual(Lines) |>
      draw
save("plots/accuracy_ancgen.png", plt)

# Comparisons across tools
accuracy_default = accuracy_default[accuracy_default.ancgen .== default_ancgen, :]
plt = data(accuracy_default) *
      mapping(:admgen,
              :mean => "Avg. accuracy",
              col = :markerset,
              color = :software) *
      visual(Lines) |>
      draw
save("plots/accuracy_software.png", plt)
