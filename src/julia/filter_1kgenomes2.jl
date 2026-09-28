

using CSV
using DataFrames
using StatsBase

# Read admixture ancestries, individuals, and global ancestries
egas = CSV.read("data/1kgenomes_AFR_EAS_EUR_chr22_maf005.3.Q", DataFrame, header = ["pop1", "pop2", "pop3"])
individuals = CSV.read("data/1kgenomes_AFR_EAS_EUR_chr22_maf005.fam", DataFrame, header = ["individual", "omit1", "omit2", "omit3", "omit4", "omit5"])[:, [:individual]]
tgas = CSV.read("data/1kgenomes_all_individuals.tsv", DataFrame, header = ["individual", "population", "superpopulation", "omit6", "omit7"])[:, 1:3]

# Combine individuals with their ancestry estimates
individuals.pop1 = egas.pop1
individuals.pop2 = egas.pop2
individuals.pop3 = egas.pop3

# Combine individuals with declared global ancestries
gas = leftjoin(individuals, tgas, on = :individual)

# Determine population numbering
sgas = combine(groupby(gas, :superpopulation), [:pop1, :pop2, :pop3] .=> mean .=> [:pop1, :pop2, :pop3])
popmatch = sort(combine(groupby(stack(sgas, [:pop1, :pop2, :pop3]), :superpopulation), [:variable, :value] => ((x,y) -> x[y .== maximum(y)]) => :i), :i)
popdict = Dict(popmatch.superpopulation .=> popmatch.i)

# Assign 
gas.pop .= 0.0
for i in axes(gas, 1)
    gas.pop[i] = gas[i, popdict[gas.superpopulation[i]]]
end

# 
gas.keep = gas.pop .>= 0.99

# Print summaries
sort(combine(groupby(gas, [:superpopulation, :population]), :keep .=> [length, mean, sum]), [:superpopulation, :keep_mean])

# Print remaining
deleteat!(gas, .~gas.keep)
CSV.write("data/1kgenomes_AFR_EAS_EUR_individuals_ga99p.txt", gas[:,[:individual]], header = false)


