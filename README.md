# Comparison between LocalAncestry.jl and existing local-ancestry estimating softwares

## Introduction

## How to run

### Install pixi environment
```{bash}
pixi install
```

### Run snakemake
```{bash}
mkdir -p logs/snakemake
nohup pixi run snakemake --profile profiles/ghpcaudk  > logs/snakemake/$(date +%Y%m%d_%H%M%S).txt 2>&1 &
```


# Set tmp directory!