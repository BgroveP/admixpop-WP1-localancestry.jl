
# Loac packages
using CSV
using DataFrames
using StatsBase

# Input
defaultt = 0.6

# 
softwares = readdir("benchmarks")
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
CSV.write("benchmarks/summary_RAM.csv", unstacked_summary_RAM)
CSV.write("benchmarks/summary_time.csv", unstacked_summary_time)
