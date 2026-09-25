
using CSV
using DataFrames
using Base.Threads
using StatsBase
using QGIO

# Input
maf = 0.005
nreplicates = 30
vcfpaths = ["ancestral", "admixed"]
predmaps = ["hapmap3", "1kgenomes"]
indsizes = 1000*collect(2:4:10)

# Read homemade functions
include("src/DUMP.jl")

# Cartesian product for analyses
writes = expand_grid(replicate = 
                    1:nreplicates, 
                    group=vcfpaths, 
                    map=predmaps, 
                    nancestral=3*indsizes, 
                    nadmixed=indsizes)

# Produce genomic maps
## Get founder allele frequencies
founderHT, _ = QGIO.haplotypes("data/extracted/_all.vcf.gz")
founderp = mean.(eachcol(founderHT))
genomic_map = "$(public_home)/data/human/genomicmap/1kgenomes_chr22_all.csv"
gmap = CSV.read(genomic_map, DataFrame)
gmap[:, "p"] = min.(founderp, 1 .- founderp)
select!(gmap, Not(["reference", "alternative"]))
gmap[:, "1kgenomes"] .= true

## Subset mark hapmap3 loci
gmaphapmap3 = CSV.read("$(public_home)/data/human/genomicmap/hapmap3_chr22_all.csv", DataFrame)
select!(gmaphapmap3, Not(["reference", "alternative", "position(cM)", "rate(cM/Mb)"]))
gmaphapmap3.hapmap3 .= true
leftjoin!(gmap, unique(gmaphapmap3), on=names(gmap)[1:2])
gmap.hapmap3[ismissing.(gmap.hapmap3)] .= false

# Write maps
gmap.mafok = gmap.p .>= maf
deleteat!(gmap, .~gmap.mafok)
gmap.id = 1:nrow(gmap)
gmap.hapmap3 = Bool.(gmap.hapmap3)
gmap.hapmapid = cumsum(gmap.hapmap3)
CSV.write("$(workdir)/map_1kgenomes_flare.ssv", gmap[:, ["chromosome", "id", "position(cM)", "position(bp)"]], delim=" ", header=false)
CSV.write("$(workdir)/map_1kgenomes_rfmix.ssv", gmap[:, ["chromosome", "position(bp)", "position(cM)"]], delim='\t', header=false)
CSV.write("$(workdir)/map_1kgenomes_recombmix.ssv", gmap[:, ["chromosome", "position(bp)", "p", "position(cM)"]], delim='\t', header=false)
CSV.write("$(workdir)/map_hapmap3_flare.ssv", gmap[gmap.hapmap3, ["chromosome", "hapmapid", "position(cM)", "position(bp)"]], delim=" ", header=false)
CSV.write("$(workdir)/map_hapmap3_rfmix.ssv", gmap[gmap.hapmap3, ["chromosome", "position(bp)", "position(cM)"]], delim='\t', header=false)
CSV.write("$(workdir)/map_hapmap3_recombmix.ssv", gmap[gmap.hapmap3, ["chromosome", "position(bp)", "p", "position(cM)"]], delim='\t', header=false)

#
i = 1
allloci = QGIO.loci("$(public_home)/simulated/human/vcf/$(scenarios.name[i])_$(scenarios.chromosome[i])_founder_all.vcf.gz")[:, 1:3]
allloci[:, "position(bp)"] = allloci[:, "position"]

tmploci = leftjoin(gmap, allloci, on=["chromosome", "position(bp)"])
tmploci.identifier = string.(tmploci.identifier)
tmploci.position = Int.(tmploci.position)

# Write all files
for i in axes(scenarios, 1)

    # In Paths
    target_path = "$(public_home)/simulated/human/vcf/$(scenarios.name[i])_$(scenarios.chromosome[i])_admixed_all.vcf.gz"
    current_path = "$(public_home)/simulated/human/vcf/$(scenarios.name[i])_$(scenarios.chromosome[i])_ancestral_all.vcf.gz"
    founder_path = "$(public_home)/simulated/human/vcf/$(scenarios.name[i])_$(scenarios.chromosome[i])_founder_all.vcf.gz"
    origins_path = "$(public_home)/simulated/human/global ancestries/$(scenarios.name[i])_all.csv"
    true_origins = "$(public_home)/simulated/human/localancestry/$(scenarios.name[i])_admixed.csv"

    # Read
    ga = CSV.read(origins_path, DataFrame, types=String)
    CSV.write("$(workdir)/ga_rfmix_$(i).tsv", ga, delim='\t', header=false)
    CSV.write("$(workdir)/ga_flare_$(i).ssv", ga, delim=" ", header=false)
    CSV.write("$(workdir)/ga_own_$(i).csv", ga)

    for j in axes(writes, 1)
        println("printing:")
        println(writes[j,:])
        if writes.group[j] != "founder"
            tga = vcf_get_n_youngest(ga, writes.indsize[j], writes.group[j])
            vcfpath = writes.group[j] == "current" ? current_path : target_path
            QGIO.subset(vcfpath, "$(workdir)/$(writes.group[j])_$(writes.m[j])_$(writes.indsize[j])_$(i).vcf.gz", samples=tga.individual, loci=tmploci[tmploci[!, writes.m[j]], :])

        else
            QGIO.subset(founder_path, "$(workdir)/$(writes.group[j])_$(writes.m[j])_all_$(i).vcf.gz", loci=tmploci[tmploci[!, writes.m[j]], :])
        end
    end
end
