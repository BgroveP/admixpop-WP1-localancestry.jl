
# Agreement between LocalAncestry.jl and the published 1000 Genomes Phase 3 local ancestry calls
# for one chromosome, evaluated at the loci of the target VCF.
#
# Usage: julia --project=. src/julia/realdata_concordance.jl <est.csv> <target.vcf.gz> <bed dir> <samples.txt> <chromosome> <out.csv>
#
# Published calls: <bed dir>/<sample>_A.bed and _B.bed, tab-separated, no header:
#   chromosome, start (bp, 0-based), end (bp), ancestry, start (cM), end (cM)
# Ancestries other than AFR/EUR (UNK = posterior < 0.9, NAT, centromere, miscall) are excluded.
#
# Output: one row per individual with counts, summed later by realdata_summary.jl:
#   hap_*  haplotype level: LocalAncestry haplotype 1/2 vs published _A/_B. The published
#          (PopPhased) calls were phase-corrected, so haplotype order may differ from the VCF.
#   dip_*  diploid level: number of AFR copies at a locus (0, 1, 2); independent of phasing.

using CSV
using DataFrames
using QGIO

estpath, vcfpath, beddir, samplespath, chrom, outpath = ARGS

const NONE = Int8(0)   # no call
const AFR = Int8(1)
const EUR = Int8(2)
const OTHER = Int8(3)  # UNK, NAT, centromere, miscall, ...
code(x) = ismissing(x) || isempty(x) ? NONE : x == "AFR" ? AFR : x == "EUR" ? EUR : OTHER
normchrom(c) = replace(string(c), r"^chr" => "")
valid(x) = (x == AFR) | (x == EUR)

# Evaluation loci
pos = sort(unique(QGIO.loci(vcfpath).position))
L = length(pos)
println("Chromosome $(chrom): $(L) loci")

# Label the loci within [lo, hi] (bp, inclusive)
function paint!(v, pos, lo, hi, val)
    a = searchsortedfirst(pos, lo)
    b = searchsortedlast(pos, hi)
    a <= b && (v[a:b] .= val)
    return v
end

# LocalAncestry.jl estimates
est = CSV.read(estpath, DataFrame, types = Dict(:individual => String, :chromosome => String, :basepairs => String, :ancestry => String))
est.lo = parse.(Int, replace.(est.basepairs, r"\:.+" => ""))
est.hi = parse.(Int, replace.(est.basepairs, r".+\:" => ""))
estgroups = groupby(est, [:individual, :haplotype])

function estimated(ind, hap)
    v = zeros(Int8, L)
    key = (individual = ind, haplotype = hap)
    haskey(estgroups, key) || (@warn "No LocalAncestry.jl estimate for $(ind) haplotype $(hap)"; return v)
    for r in eachrow(estgroups[key])
        paint!(v, pos, r.lo, r.hi, code(r.ancestry))
    end
    return v
end

# Published calls
function published(ind, suffix)
    v = zeros(Int8, L)
    path = joinpath(beddir, "$(ind)_$(suffix).bed")
    bed = CSV.read(path, DataFrame, header = false, delim = '\t', types = Dict(1 => String, 4 => String))
    for r in eachrow(bed)
        normchrom(r[1]) == chrom || continue
        paint!(v, pos, r[2] + 1, r[3], code(r[4]))   # 0-based start -> 1-based
    end
    return v
end

# Count agreement for one individual (t: published A/B, e: estimated 1/2)
function counts(t, e)
    hap_truth_valid = hap_truth_excluded = hap_assigned = hap_correct = 0
    dip_valid = dip_correct = dip_abserr = afr_true = afr_est = 0

    # Haplotype level
    for h in 1:2, k in eachindex(t[h])
        th, eh = t[h][k], e[h][k]
        if valid(th)
            hap_truth_valid += 1
            if valid(eh)
                hap_assigned += 1
                hap_correct += th == eh
            end
        else
            hap_truth_excluded += 1
        end
    end

    # Diploid level (AFR dosage)
    for k in eachindex(t[1])
        (valid(t[1][k]) & valid(t[2][k]) & valid(e[1][k]) & valid(e[2][k])) || continue
        dt = (t[1][k] == AFR) + (t[2][k] == AFR)
        de = (e[1][k] == AFR) + (e[2][k] == AFR)
        dip_valid += 1
        dip_correct += dt == de
        dip_abserr += abs(dt - de)
        afr_true += dt
        afr_est += de
    end
    return (; hap_truth_valid, hap_truth_excluded, hap_assigned, hap_correct,
              dip_valid, dip_correct, dip_abserr, afr_true, afr_est)
end

samples = String.(filter(!isempty, strip.(readlines(samplespath))))
rows = DataFrame([(individual = ind, chromosome = chrom, nloci = L,
                   counts((published(ind, "A"), published(ind, "B")), (estimated(ind, 1), estimated(ind, 2)))...)
                  for ind in samples])

CSV.write(outpath, rows)
println("Diploid concordance chromosome $(chrom): ", sum(rows.dip_correct) / max(1, sum(rows.dip_valid)))
