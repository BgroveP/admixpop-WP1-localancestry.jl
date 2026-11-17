
# Summarise per-locus accuracies across replicates.
# Usage: julia -t <threads> --project=. src/julia/result_accuracies.jl <tool1,tool2,...> <default threshold>

using Base.Threads
using CSV
using DataFrames
using Statistics
using StatsBase

# Input
softwares = length(ARGS) >= 1 ? String.(split(ARGS[1], ",")) :
            ["flare", "recombmix", "gnomix", "localancestryjl", "rfmix", "sparsepainter"]
defaultt = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.66
softwares = softwares[isdir.("results/accuracy/" .* softwares)]
nsoftwares = length(softwares)
mkpath("results/summary")

# File names: <replicate>_<markerset>_<ancgen>_<admgen>[_<threshold>].csv
const NAME_RE = r"^(\d+)_([A-Za-z]+)_(\d+)_(\d+)(?:_([0-9.]+?))?\.csv$"

# Read accuracies
o = Vector{Vector{DataFrame}}(undef, nsoftwares)
@threads for i in axes(softwares, 1)
    s = softwares[i]
    files = filter(f -> occursin(NAME_RE, basename(f)), readdir("results/accuracy/$(s)", join = true))
    o[i] = Vector{DataFrame}(undef, length(files))

    for (ji, j) in enumerate(files)
        df = CSV.read(j,
                      DataFrame,
                      types = [String, Int, Float64, Bool, Float64, Float64, Float64],
                      stripwhitespace = true,
                      strict = true)
        mt = match(NAME_RE, basename(j))

        df.software  .= s
        df.replicate .= parse(Int, mt[1])
        df.markerset .= mt[2]
        df.ancgen    .= parse(Int, mt[3])
        df.admgen    .= parse(Int, mt[4])
        # Only LocalAncestry.jl has a threshold; other tools get the default so tables align
        df.threshold .= isnothing(mt[5]) ? defaultt : parse(Float64, mt[5])
        o[i][ji] = df
    end
end

# Collect into one dataframe
olong = reduce(vcat, reduce(vcat, o))

# Omit locus that is never estimated (50799821) and the one never estimated for rfmix and gnomix (50780578)
olong.keep = (olong[!, "position(bp)"] .!= 50799821) .& (olong[!, "position(bp)"] .!= 50780578)
olong.keep = olong.keep .& ((olong.array .& (olong.markerset .== "array")) .| (olong.markerset .== "sequence"))
deleteat!(olong, .~olong.keep)

# Which analyses are missing (scenario x replicate combinations not seen for a tool)
oamount = unique(olong[:, ["software", "ancgen", "admgen", "markerset", "threshold", "replicate"]])
oamount.n .= true
oamountseen = stack(unstack(oamount, :software, :n), softwares, variable_name = :software, value_name = :seen)
deleteat!(oamountseen, .~ismissing.(oamountseen.seen))
deleteat!(oamountseen, (oamountseen.threshold .!= defaultt) .& (oamountseen.software .!= "localancestryjl"))
CSV.write("results/summary/missing_accuracy.csv", oamountseen)
nrow(oamountseen) > 0 && @warn "$(nrow(oamountseen)) tool/scenario/replicate combinations are missing; see results/summary/missing_accuracy.csv"

# Mean accuracies
olong_means = combine(groupby(olong, ["software", "ancgen", "admgen", "markerset", "array", "threshold"]), :accuracy => mean => :accuracy)
owide_means = sort(unstack(olong_means, :software, :accuracy), [:threshold, :array, :ancgen, :admgen, :markerset])
CSV.write("results/summary/accuracy.csv", owide_means)

# Standard deviation across replicates of the replicate mean accuracy
olong_means2 = combine(groupby(olong, ["software", "ancgen", "admgen", "markerset", "array", "threshold", "replicate"]), :accuracy => mean => :accuracy)
olong_sds = combine(groupby(olong_means2, ["software", "ancgen", "admgen", "markerset", "array", "threshold"]), :accuracy => std => :sd)
owide_sds = sort(unstack(olong_sds, :software, :sd), [:threshold, :array, :ancgen, :admgen, :markerset])
CSV.write("results/summary/accuracy_sds.csv", owide_sds)
