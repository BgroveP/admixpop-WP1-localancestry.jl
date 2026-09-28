
using CSV, DataFrames, LocalAncestry

packageroot = dirname(dirname(pathof(LocalAncestry)))
x = localancestry("$(packageroot)/data/reference.vcf.gz", 
                  "$(packageroot)/data/target.vcf.gz", 
                  "$(packageroot)/data/globalancestries.csv")

