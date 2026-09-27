
# Load the Pkg manager
using Pkg

# Install official packages
Pkg.add(["DataFrames", "CSV", "StatsBase"])
Pkg.add(name="LocalAncestry", version="0.3.0")

# Install a VCF parser from github
Pkg.add(url= "https://github.com/BgroveP/QGIO.jl", rev="c196f98fbffc372913de72b4c22af4207cba2f0f")

