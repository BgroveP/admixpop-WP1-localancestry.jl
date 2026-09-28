
using CSV
using DataFrames
using StatsBase
using QGIO

# Input
maf = 0.005
nreplicates = 30
vcfpaths = ["ancestral", "admixed"]
predmaps = ["array", "sequence"]
generations =  ["6", "8", "10"] .* "_10"

# Produce genomic maps
## Get founder allele frequencies for sequence loci
founderHT, _ = QGIO.haplotypes("data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz")
founderp = mean.(eachcol(founderHT))
genomic_map = QGIO.loci("data/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz")
genomic_map[:, "p"] = min.(founderp, 1 .- founderp)
genomic_map[:, "sequence"] = genomic_map[:, "p"] .>= maf

## Join with hapmap3 loci
genomic_map_hapmap3 = CSV.read("data/hapmap3.hg38.bed", DataFrame, header = [:chromosome, :position, :omit1])
genomic_map_hapmap3.array .= true
leftjoin!(genomic_map, unique(genomic_map_hapmap3), on= [:chromosome, :position])
genomic_map.array[ismissing.(genomic_map.array)] .= false
genomic_map.array = genomic_map.array .& genomic_map.sequence

# Remove loci with maf lower than threshold
deleteat!(genomic_map, genomic_map.p .< maf)

# Get locusid for sequence and array marker sets
genomic_map.sequenceid = 1:nrow(genomic_map)
genomic_map.arrayid = cumsum(genomic_map.array)

# Insert recombination rates used for simulation
genomic_map.distance =  genomic_map.position * 10^-6
genomic_map.dummy = genomic_map.position .-1

# Write maps for each combination between software and marker set
for s in predmaps
    CSV.write("data/flare_$(s).ssv", genomic_map[genomic_map[:, s], ["chromosome", "$(s)id", "distance", "position"]], delim=" ", header=false)
    CSV.write("data/rfmix_$(s).tsv", genomic_map[genomic_map[:, s], ["chromosome", "position", "distance"]], delim='\t', header=false)
    CSV.write("data/recombmix_$(s).tsv", genomic_map[genomic_map[:, s], ["chromosome", "position", "p", "distance"]], delim='\t', header=false)
    CSV.write("data/bcftools_$(s).tsv", genomic_map[genomic_map[:, s], ["chromosome", "position"]], delim='\t', header=false)
end
