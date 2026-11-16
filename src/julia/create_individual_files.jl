
using CSV
using DataFrames

datapath = "data"
replicatepaths = readdir(datapath, join = true)
replicatepaths = replicatepaths[contains.(replicatepaths, "replicate_" )]
# Generations (first generation of each panel), e.g. "6,8,10"
steps = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], ",")) : [6, 8, 10]

for r in replicatepaths

    globalancestries = CSV.read("$(r)/global_ancestries.csv", DataFrame)
    deleteat!(globalancestries, contains.(globalancestries.population, "anc", ))
    transform!(groupby(globalancestries, :population), :population => (x -> 5 .+ Int.(ceil.(collect(1:length(x)) ./ 2000))) => :generation)

    # Write global ancestries
    CSV.write("$(r)/global_ancestries_rfmix.tsv", globalancestries[globalancestries.population .!= "admixed", [:individual, :population]], delim='\t', header=false)
    CSV.write("$(r)/global_ancestries_flare.ssv", globalancestries[globalancestries.population .!= "admixed", [:individual, :population]], delim=" ", header=false)

    for s in steps
        tmpga = copy(globalancestries)
        deleteat!(tmpga, tmpga.generation .< s)
        # ancestral
        CSV.write("$(r)/individuals_ancestral_gen$(s)_10.txt", select(tmpga, :individual)[tmpga.population .!= "admixed",:], writeheader = false)
        # admixed
        CSV.write("$(r)/individuals_admixed_gen$(s)_10.txt", select(tmpga, :individual)[tmpga.population .== "admixed",:], writeheader = false)
        # Write global ancestries for gnomix
        CSV.write("$(r)/global_ancestries_gnomix_$(s).tsv", tmpga[tmpga.population .!= "admixed", [:individual, :population]], delim='\t', header=false)

    end
end
