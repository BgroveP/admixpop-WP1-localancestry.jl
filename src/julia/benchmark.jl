
# Loac packages
using CSV
using DataFrames
using StatsBase

# Input
defaultt = 0.6

# 
softwares = readdir("benchmarks")
softwares = softwares[.~contains.(softwares, ".csv")]
nsoftwares = length(softwares)

# Read
o = Vector{DataFrame}(undef, nsoftwares)
for (i,s) in enumerate(softwares)
    b = readdir("benchmarks/$(s)", join = true)
    o[i] = CSV.read(b, DataFrame)
    o[i].file = replace.(b, "benchmarks/" =>"")
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

# Summary
# Number of finishes analyses 
n = sort(unstack(combine(groupby(o, [:software, :ancgen, :admgen, :markerset, :threshold]), :s => length => :n), :software, :n), [:markerset, :ancgen, :admgen])
deleteat!(n, n.threshold .!= defaultt)

# Statistics
summary = sort(combine(groupby(o, [:software, :ancgen, :admgen, :markerset, :threshold]), [:s, :max_rss] .=> mean), :s_mean)
deleteat!(summary, summary.threshold .!= defaultt)
unstacked_summary_time = unstack(summary[:, Not(:max_rss_mean)], :software, :s_mean)
unstacked_summary_RAM = unstack(summary[:, Not(:s_mean)], :software, :max_rss_mean)
CSV.write("benchmarks/summary_RAM_mean.csv", unstacked_summary_RAM)
CSV.write("benchmarks/summary_time_mean.csv", unstacked_summary_time)

# Maxes for slurm resource allocations
summary = sort(combine(groupby(o, [:software, :ancgen, :admgen, :markerset, :threshold]), [:s, :max_rss] .=> maximum .=> [:s, :max_rss]), :s)
deleteat!(summary, summary.threshold .!= defaultt)
unstacked_summary_time = unstack(summary[:, Not(:max_rss)], :software, :s)
unstacked_summary_RAM = unstack(summary[:, Not(:s)], :software, :max_rss)
unstacked_summary_time[:, 5:end] = unstacked_summary_time[:, 5:end] ./ 60
unstacked_summary_RAM[:, 5:end] = unstacked_summary_RAM[:, 5:end] ./ 1000
CSV.write("benchmarks/summary_RAM_max.csv", unstacked_summary_RAM)
CSV.write("benchmarks/summary_time_max.csv", unstacked_summary_time)

