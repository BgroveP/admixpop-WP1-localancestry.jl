
# Summarise runtime and memory from Snakemake benchmark files.
# Usage: julia --project=. src/julia/result_benchmark.jl <tool1,tool2,...> <default threshold>

using CSV
using DataFrames
using StatsBase

# Input
softwares = length(ARGS) >= 1 ? String.(split(ARGS[1], ",")) : filter(isdir, readdir("benchmarks", join = true)) .|> basename
defaultt = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.66
softwares = softwares[isdir.("benchmarks/" .* softwares)]
mkpath("results/summary")

# File names: <replicate>_<markerset>_<ancgen>_<admgen>[_<threshold>].tsv
const NAME_RE = r"^(\d+)_([A-Za-z]+)_(\d+)_(\d+)(?:_([0-9.]+?))?\.tsv$"

# Read (only the columns that are used; `h:m:s` can exceed 24 h and is not parsed)
o = DataFrame[]
for s in softwares
    files = filter(f -> occursin(NAME_RE, basename(f)), readdir("benchmarks/$(s)", join = true))
    for f in files
        df = CSV.read(f, DataFrame, select = ["s", "max_rss"], types = Float64)
        mt = match(NAME_RE, basename(f))
        df.software  .= s
        df.replicate .= parse(Int, mt[1])
        df.markerset .= mt[2]
        df.ancgen    .= parse(Int, mt[3])
        df.admgen    .= parse(Int, mt[4])
        df.threshold .= isnothing(mt[5]) ? defaultt : parse(Float64, mt[5])
        push!(o, df)
    end
end
o = reduce(vcat, o)
deleteat!(o, o.threshold .!= defaultt)

# Number of finished runs per scenario
println(sort(unstack(combine(groupby(o, [:software, :ancgen, :admgen, :markerset]), :s => length => :n), :software, :n), [:markerset, :ancgen, :admgen]))

# Statistics
sscale = 1 / 60       # seconds -> minutes
memscale = 1 / 1024   # MB -> GB
ndigits = 2
groupcols = [:software, :ancgen, :admgen, :markerset, :threshold]

function summarise(f)
    out = sort(combine(groupby(o, groupcols), [:s, :max_rss] .=> f .=> [:s, :mem]), :s)
    out.s = round.(out.s * sscale, digits = ndigits)
    out.mem = round.(out.mem * memscale, digits = ndigits)
    return out
end

means = summarise(mean)
CSV.write("results/summary/time_use_mean.csv", unstack(means[:, Not(:mem)], :software, :s))
CSV.write("results/summary/RAM_use_mean.csv", unstack(means[:, Not(:s)], :software, :mem))

sds = summarise(std)
CSV.write("results/summary/time_use_sds.csv", unstack(sds[:, Not(:mem)], :software, :s))
CSV.write("results/summary/RAM_use_sds.csv", unstack(sds[:, Not(:s)], :software, :mem))

# LaTeX-ready table (mean_{sd}) for the largest reference panel
tbl_ancgen = maximum(o.ancgen)
meansandsds = leftjoin(means, sds, on = groupcols, makeunique = true)
deleteat!(meansandsds, meansandsds.ancgen .!= tbl_ancgen)
meansandsds.time = "\$" .* string.(meansandsds.s) .* "_{" .* string.(meansandsds.s_1) .* "}\$"
meansandsds.RAM = "\$" .* string.(meansandsds.mem) .* "_{" .* string.(meansandsds.mem_1) .* "}\$"
meansandsdslong = stack(meansandsds, [:time, :RAM])
meansandsdslong.col = string.(meansandsdslong.admgen) .* string.(meansandsdslong.variable)

meansandsdswider = unstack(meansandsdslong[:, [:software, :markerset, :value, :col]], :col, :value)
sort!(meansandsdswider, [:markerset, :software])
admgens = sort(unique(meansandsdslong.admgen), rev = true)
colorder = [string.(admgens) .* "time"; string.(admgens) .* "RAM"]
meansandsdswider = meansandsdswider[:, [["software", "markerset"]; colorder]]
meansandsdswider.outstring = [join(coalesce.(collect(r[colorder]), "-"), " & ") for r in eachrow(meansandsdswider)]

CSV.write("results/summary/resources_table_latex.csv", meansandsdswider[:, [:software, :markerset, :outstring]])
