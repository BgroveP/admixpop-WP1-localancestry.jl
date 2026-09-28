


using CSV
using DataFrames
using QGIO

# 1000 Genomes
kgenomes = CSV.read("data/1kgenomes_all_individuals.tsv", DataFrame, delim = "\t")
deleteat!(kgenomes, .~((kgenomes.super_pop .== "EUR") .| (kgenomes.super_pop .== "AFR") .| (kgenomes.super_pop .== "EAS")))
rename!(kgenomes, [:individual, :population, :superpopulation, :sex, :rubbish1, :rubbish2])

# ISGR
isgr = CSV.read("data/all_individuals.tsv", DataFrame, delim = "\t")
deleteat!(isgr, ismissing.(isgr.superpopulation))

# VCF 
vcf = DataFrame(individual = QGIO.samples("data/1kgenomes_all_chr22.vcf.gz"))

# Intersect
out = innerjoin(innerjoin(isgr, kgenomes, on = [:individual, :superpopulation]),vcf, on = :individual)

# Write individuals file
CSV.write("data/1kgenomes_AFR_EAS_EUR_individuals.txt", out[:,[:individual]], writeheader = false)