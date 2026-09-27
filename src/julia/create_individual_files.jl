
using CSV
using DataFrames

datapath = "data/extracted"
replicatepaths = readdir(datapath, join = true)
steps = [6,8,10]

for r in replicatepaths

    globalancestries = CSV.read("$(r)/global_ancestries.csv", DataFrame)
    deleteat!(globalancestries, contains.(globalancestries.population, "anc", ))
    transform!(groupby(globalancestries, :population), :population => (x -> 5 .+ Int.(ceil.(collect(1:length(x)) ./ 2000))) => :generation)

    for s in steps
        tmpga = copy(globalancestries)
        deleteat!(tmpga, tmpga.generation .< s)
    # ancestral
        CSV.write("$(r)/individuals_ancestral_gen$(s)_10.txt", select(tmpga, :individual)[tmpga.population .!= "admixed",:], writeheader = false)
    # admixed
        CSV.write("$(r)/individuals_admixed_gen$(s)_10.txt", select(tmpga, :individual)[tmpga.population .== "admixed",:], writeheader = false)

    end
end
