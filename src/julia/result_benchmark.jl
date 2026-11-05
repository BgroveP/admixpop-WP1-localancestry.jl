
# Load packages
using CSV
using DataFrames
using StatsBase
using Dates 

# Input
defaultt = 0.6

# 
softwares = readdir("benchmarks")
softwares = softwares[.~contains.(softwares, ".csv")]
nsoftwares = length(softwares)

# Read
o = Vector{DataFrame}(undef, nsoftwares)
for (i,s) in enumerate(softwares)
    println(s)
    b = readdir("benchmarks/$(s)", join = true)
    o[i] = CSV.read(b, DataFrame, types = [Float64, Dates.Time, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64])
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
# Number of missing analyses 
n = sort(unstack(combine(groupby(o, [:software, :ancgen, :admgen, :markerset, :threshold]), :s => length => :n), :software, :n), [:markerset, :ancgen, :admgen])
deleteat!(n, n.threshold .!= defaultt)
nmissing = transform(n, [:flare, :gnomix, :localancestryjl, :recombmix, :rfmix, :sparsepainter] .=> (x -> 30 .- x) .=> [:flare, :gnomix, :localancestryjl, :recombmix, :rfmix, :sparsepainter])
Dict([:flare, :flare_clust, :gnomix, :localancestryjl, :recombmix, :rfmix, :sparsepainter] .=> sum.(eachcol(nmissing[:,5:end])))

# Statistics
sscale = 1/60 # Seconds to minutes
sdigits = 2
memscale = 1/1024 # Megabytes to Gigabytes
memdigits = 2
means = sort(combine(groupby(o, [:software, :ancgen, :admgen, :markerset, :threshold]), [:s, :max_rss] .=> mean .=> [:s, :mem]), :s)
deleteat!(means, means.threshold .!= defaultt)
means.s = round.(means.s * sscale, digits = sdigits)
means.mem = round.(means.mem * memscale, digits = memdigits)
unstacked_means_time = unstack(means[:, Not(:mem)], :software, :s)
unstacked_means_RAM = unstack(means[:, Not(:s)], :software, :mem)
CSV.write("results/summary/RAM_use_mean.csv",  unstacked_means_RAM)
CSV.write("results/summary/time_use_mean.csv", unstacked_means_time)

sds = sort(combine(groupby(o, [:software, :ancgen, :admgen, :markerset, :threshold]), [:s, :max_rss] .=> (x -> sqrt(var(x))) .=> [:s, :mem]), :s)
deleteat!(sds, sds.threshold .!= defaultt)
sds.s = round.(sds.s * sscale, digits = sdigits)
sds.mem = round.(sds.mem * memscale, digits = memdigits)
unstacked_sds_time = unstack(sds[:, Not(:mem)], :software, :s)
unstacked_sds_RAM = unstack(sds[:, Not(:s)], :software, :mem)
CSV.write("results/summary/RAM_use_sds.csv",  unstacked_sds_RAM)
CSV.write("results/summary/time_use_sds.csv", unstacked_sds_time)

# Write table
tbl_ancgen = 10
meansandsds = leftjoin(means, sds, on = names(means)[1:5], makeunique = true)
deleteat!(meansandsds, meansandsds.ancgen .!= tbl_ancgen)
meansandsds.time = "\$" .*string.(meansandsds.s) .* "_{" .* string.(meansandsds.s_1) .* "}\$"
meansandsds.RAM = "\$" .*string.(meansandsds.mem) .* "_{" .* string.(meansandsds.mem_1) .* "}\$"
meansandsdslong = stack(meansandsds, [:time, :RAM])
meansandsdslong.col = string.(meansandsdslong.admgen) .* string.(meansandsdslong.variable)

meansandsdswider = unstack(meansandsdslong[:, [:software, :markerset, :value,:col]], :col, :value)
sort!(meansandsdswider, [:markerset, :software])
colorder = ["10time", "8time", "6time", "10RAM", "8RAM", "6RAM"]
meansandsdswider = meansandsdswider[:, [["software", "markerset"]; colorder]]
meansandsdswider.outstring .= ""
for i in eachrow(meansandsdswider)
    i[end] = join(i[3:(end-1)], " & ")
end

CSV.write("tbloutstring.csv", meansandsdswider[:, [:software, :markerset, :outstring]])
