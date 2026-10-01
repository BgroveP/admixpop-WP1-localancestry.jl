
# Loac packages
using CSV
using DataFrames
using StatsBase

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

# Summary
s1 = sort(combine(groupby(o, [:software, :ancgen, :admgen, :markerset]), [:s, :max_rss] .=> mean), :s_mean)
s2 = transform(groupby(s1, [:ancgen, :admgen, :markerset]), 
[:software, :s_mean] => ((x,y) -> y ./ y[x .== "localancestryjl"]) => :s,
[:software, :max_rss_mean] => ((x,y) -> y ./ y[x .== "localancestryjl"]) => :rss)
s3 = combine(groupby(s2, :software), :s .=> [mean, median, minimum, maximum], :rss .=> [mean, median, minimum, maximum], )
