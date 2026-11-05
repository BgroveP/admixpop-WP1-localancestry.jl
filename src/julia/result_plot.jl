
using CSV
using DataFrames
using Statistics
using StatsBase
using Printf
using AlgebraOfGraphics
using CairoMakie

# Input
softwares = [:flare, :recombmix, :gnomix, :localancestryjl, :rfmix, :sparsepainter]
default_threshold = 0.6
default_ancgen = 10

# Translations
msadict = Dict(
    ["sequence", false] => "Sequence -> Sequence",
    ["sequence", true] => "Sequence -> Array",
    ["array", true] => "Array -> Array",
)

# Read the accuracies
accuracies = stack(CSV.read("results/summary/accuracy.csv", DataFrame, stripwhitespace=true, missingstring = ["NaN", ""]), softwares, variable_name=:software, value_name=:mean)
deleteat!(accuracies, ismissing.(accuracies.mean))
accuraciessd = stack(CSV.read("results/summary/accuracy_sds.csv", DataFrame, stripwhitespace=true, missingstring = ["NaN", ""]), softwares, variable_name=:software, value_name=:sd)
deleteat!(accuraciessd, ismissing.(accuraciessd.sd))
leftjoin!(accuracies, accuraciessd, on = names(accuracies)[1:(end-1)])

# Translate columns
accuracies.markerset = String.(accuracies.markerset)
for i in axes(accuracies, 1)
    accuracies.markerset[i] = get(msadict, [accuracies.markerset[i], accuracies.array[i]], accuracies.markerset[i])
end

# Little accuracy increase from larger reference panel
accuracy_default = accuracies[accuracies.threshold .== default_threshold,:]
plt = data(accuracy_default) *
          mapping(:ancgen,
                  :mean => "Avg. accuracy",
                  col = :admgen,
                  color = :software,
                  group = (:software, :markerset, :admgen)) *
          visual(Lines) |> 
          draw
save("plots/accuracy_ancgen.png", plt)

# Comparisons across softwares
accuracy_default = accuracy_default[accuracy_default.ancgen .== default_ancgen,:]
plt = data(accuracy_default) *
          mapping(:admgen ,
                  :mean => "Avg. accuracy",
                  col = :markerset,
                  color = :software) *
          visual(Lines) |> 
          draw
save("plots/accuracy_software.png", plt)

"""
    plot_accuracy_patterns(df; threshold = 0.6, path = "accuracy_patterns_0.6.png")

AlgebraOfGraphics figure: mean accuracy per method, split by marker set.
"""
function plot_accuracy_patterns(df::AbstractDataFrame; threshold = 0.6,
                                path = "accuracy_patterns_0.6.png")
    sub = filter(:threshold => ==(threshold), df)

    # long format: one row per (config, method) with non-missing accuracy
    long = stack(sub, METHODS; variable_name = :method, value_name = :accuracy) |>
           t -> dropmissing(t, :accuracy)

    summary = groupby(long, [:markerset, :method]) |>
              g -> combine(g, :accuracy => mean => :mean_accuracy)

    plt = data(summary) *
          mapping(:method => "Method",
                  :mean_accuracy => "Mean accuracy",
                  color = :markerset,
                  dodge = :markerset) *
          visual(BarPlot; dodge_gap = 0.1)

    fg = draw(plt;
              axis = (title = "Accuracy by method at threshold = $threshold",
                      limits = (nothing, (0.80, 0.92))))
    savefig(fg, path)
    println("figure written to $path")
    return fg
end

# --- Run ------------------------------------------------------------------
function main()
    df = CSV.read(IOBuffer(BENCH_CSV), DataFrame; missingstring = [""])
    accuracy_patterns(df; threshold = 0.6)
    plot_accuracy_patterns(df; threshold = 0.6)
end

abspath(PROGRAM_FILE) == @__FILE__ && main()