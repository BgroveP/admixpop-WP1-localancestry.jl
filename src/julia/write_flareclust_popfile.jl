

# Import libraries
using CSV
using DataFrames
using QGIO

# Get commandline arguments
vcfpath = ARGS[1]
poppath = ARGS[2]
outpath = ARGS[3]

# Read data from files
poplabels = CSV.read(poppath, DataFrame)
samples = DataFrame(individual = QGIO.samples(vcfpath))

# Combine information to only have information of individuals in the vcf file.
leftjoin!(samples, poplabels, on = :individual)

# Write the resulting information
CSV.write(outpath, samples, header = false, delim = ' ')
