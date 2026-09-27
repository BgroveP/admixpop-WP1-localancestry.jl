
using CSV
using DataFrames
using Base.Threads
using StatsBase
using QGIO

# Input
maf = 0.005
nreplicates = 30
vcfpaths = ["ancestral", "admixed"]
predmaps = ["array", "sequence"]
generations =  ["6", "8", "10"] .* "_10"

#
include("DUMP.jl")
# Cartesian product for analyses
writes = expand_grid(replicate = 
                    1:nreplicates, 
                    group=vcfpaths, 
                    map=predmaps, 
                    ancgens=generations, 
                    admgens=generations)

# Produce genomic maps
## Get founder allele frequencies for sequence loci
founderHT, _ = QGIO.haplotypes("data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz")
founderp = mean.(eachcol(founderHT))
genomic_map = QGIO.loci("data/1kgenomes/1kgenomes_AFR_EAS_EUR_chr22_maf005_ga99p.vcf.gz")
genomic_map[:, "p"] = min.(founderp, 1 .- founderp)
genomic_map[:, "sequence"] = genomic_map[:, "p"] .>= maf

## Join with hapmap3 loci
genomic_map_hapmap3 = CSV.read("data/maps/hapmap3.hg38.bed", DataFrame, header = [:chromosome, :position, :omit1])
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
    CSV.write("data/maps/flare_$(s).ssv", genomic_map[genomic_map[:, s], ["chromosome", "$(s)id", "distance", "position"]], delim=" ", header=false)
    CSV.write("data/maps/rfmix_$(s).tsv", genomic_map[genomic_map[:, s], ["chromosome", "position", "distance"]], delim='\t', header=false)
    CSV.write("data/maps/recombmix_$(s).tsv", genomic_map[genomic_map[:, s], ["chromosome", "position", "p", "distance"]], delim='\t', header=false)
    CSV.write("data/maps/bcftools_$(s).bed", genomic_map[genomic_map[:, s], ["chromosome", "dummy", "position"]], delim='\t', header=false)
end

# Write genomic maps for subsetting


# Write global ancestries
for i in 1:nreplicates

    # In Paths
    globalancestries_path = "data/extracted/replicate_$(i)/global_ancestries.csv"

    # Read the global ancestry file
    ga = CSV.read(globalancestries_path, DataFrame, types=String)

    # Write the global ancestry files for each software
    CSV.write("data/extracted/replicate_$(i)/rfmix_globalancestry.tsv", ga, delim='\t', header=false)
    CSV.write("data/extracted/replicate_$(i)/flare_globalancestry.ssv", ga, delim=" ", header=false)
    CSV.write("data/extracted/replicate_$(i)/localancestryjl_globalancestry.csv", ga)

end

# Subset the vcf files
writes = expand_grid(replicate = 1:nreplicates, map=predmaps, gens=generations)


@threads for (r, m, g) in eachrow(writes)

    println("rep: $(r) map: $(m) gens: $(g)")
    flush(stdout)

    # In Paths
    ancestralvcf_path = "data/extracted/replicate_$(r)/ancestral.vcf.gz"
    admixedvcf_path = "data/extracted/replicate_$(r)/admixed.vcf.gz"

    # Read lists of focal individuals
    individuals_ancestral_path = "data/extracted/replicate_$(r)/individuals_ancestral_gen$(g).txt"
    individuals_admixed_path = "data/extracted/replicate_$(r)/individuals_admixed_gen$(g).txt"
    
    # Subset vcf files
    run(`bcftools view -R data/maps/bcftools_$(m).bed -S $(individuals_ancestral_path) $(ancestralvcf_path) -Oz -o data/extracted/replicate_$(r)/ancestral_$(m)_$(g).vcf.gz`)
    run(`bcftools view -R data/maps/bcftools_$(m).bed -S $(individuals_admixed_path) $(admixedvcf_path) -Oz -o data/extracted/replicate_$(r)/admixed_$(m)_$(g).vcf.gz`)

    # Index the vcf files
    run(`tabix -p vcf data/extracted/replicate_$(r)/ancestral_$(m)_$(g).vcf.gz`)
    run(`tabix -p vcf data/extracted/replicate_$(r)/admixed_$(m)_$(g).vcf.gz`)
    
end
