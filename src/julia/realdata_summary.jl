
# Genome-wide summary of LocalAncestry.jl on real data vs. the published 1000 Genomes calls.
# Usage: julia --project=. src/julia/realdata_summary.jl <target population> <chromosomes, comma-separated>
#
# Writes results/realdata/
#   summary.csv      concordance per chromosome and genome-wide ("all")
#   individuals.csv  genome-wide AFR proportion (published vs. LocalAncestry.jl) and concordance per individual
#   resources.csv    LocalAncestry.jl runtime and peak memory per chromosome and genome-wide

using CSV
using DataFrames
using Statistics

target = ARGS[1]
chroms = String.(split(ARGS[2], ","))
outdir = "results/realdata"

ind = reduce(vcat, [CSV.read("$(outdir)/$(target)_chr$(c)_individuals.csv", DataFrame,
                             types = Dict(:individual => String, :chromosome => String)) for c in chroms])

# --- Concordance ---
function concordance(df)
    hap_loci = sum(df.hap_truth_valid) + sum(df.hap_truth_excluded)
    return (
        nindividuals = length(unique(df.individual)),
        nloci = sum(unique(df[:, [:chromosome, :nloci]]).nloci),
        # Diploid AFR dosage (phase-independent): main metric
        diploid_concordance = sum(df.dip_correct) / sum(df.dip_valid),
        diploid_mean_abs_error = sum(df.dip_abserr) / sum(df.dip_valid),
        # Haplotype level (sensitive to phase differences between the VCF and the published calls)
        haplotype_concordance = sum(df.hap_correct) / sum(df.hap_assigned),
        # Fraction of haplotype-loci with a published AFR/EUR call that LocalAncestry.jl assigned
        haplotype_coverage = sum(df.hap_assigned) / sum(df.hap_truth_valid),
        # Fraction of haplotype-loci excluded because the published call is not AFR/EUR (UNK, NAT, ...)
        published_excluded = sum(df.hap_truth_excluded) / hap_loci,
    )
end

summary = DataFrame([(chromosome = k.chromosome, concordance(g)...) for (k, g) in pairs(groupby(ind, :chromosome; sort = false))])
push!(summary, (chromosome = "all", concordance(ind)...))

# --- Per individual, genome-wide ---
perind = combine(groupby(ind, :individual),
                 [:afr_true, :dip_valid] => ((a, n) -> sum(a) / (2 * sum(n))) => :afr_published,
                 [:afr_est, :dip_valid] => ((a, n) -> sum(a) / (2 * sum(n))) => :afr_localancestryjl,
                 [:dip_correct, :dip_valid] => ((c, n) -> sum(c) / sum(n)) => :diploid_concordance)
# Global AFR proportions are only meaningful genome-wide ("all" row)
summary.global_afr_correlation = Vector{Union{Missing,Float64}}(missing, nrow(summary))
summary.global_afr_mean_abs_diff = Vector{Union{Missing,Float64}}(missing, nrow(summary))
summary[end, :global_afr_correlation] = cor(perind.afr_published, perind.afr_localancestryjl)
summary[end, :global_afr_mean_abs_diff] = mean(abs.(perind.afr_published .- perind.afr_localancestryjl))

# --- Resources (Snakemake benchmark files) ---
res = DataFrame()
for c in chroms
    b = CSV.read("benchmarks/realdata/localancestryjl/$(target)_chr$(c).tsv", DataFrame, select = ["s", "max_rss"], types = Float64)
    push!(res, (chromosome = c, minutes = b.s[1] / 60, max_rss_gb = b.max_rss[1] / 1024))
end
push!(res, (chromosome = "all", minutes = sum(res.minutes), max_rss_gb = maximum(res.max_rss_gb)))

CSV.write("$(outdir)/summary.csv", summary)
CSV.write("$(outdir)/individuals.csv", perind)
CSV.write("$(outdir)/resources.csv", res)

println(summary)
println(res)
