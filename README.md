# Comparison between LocalAncestry.jl and existing local ancestry inference tools

A Snakemake pipeline that compares the accuracy and computational cost of local ancestry
inference (LAI) tools on simulated three-way admixed data (AFR, EAS and EUR, chromosome 22):

| Tool | Config name |
|---|---|
| [LocalAncestry.jl](https://github.com/BgroveP/LocalAncestry.jl) | `localancestryjl` |
| [FLARE](https://github.com/browning-lab/flare) | `flare` (and optional `flareclust`) |
| [RFMix v2](https://github.com/slowkoni/rfmix) | `rfmix` |
| [Recomb-Mix](https://github.com/ucfcbb/Recomb-Mix) | `recombmix` |
| [SparsePainter](https://github.com/YaolingYang/SparsePainter) | `sparsepainter` |
| [G-Nomix](https://github.com/AI-sandbox/gnomix) | `gnomix` |

Each tool is run for every combination of
- **replicate** (`replicates`, simulated datasets from Zenodo),
- **marker set** (`markersets`: `array` = HapMap3 positions, `sequence` = all loci with MAF ≥ 0.005),
- **reference panel** (`generations` → `g1`: individuals from generation g1–10), and
- **target panel** (`generations` → `g2`: admixed individuals from generation g2–10).

## Requirements

- Linux x86-64 (the pixi environment is defined for `linux-64` only)
- [pixi](https://pixi.sh)
- `git` and `curl` on the PATH
- Lots of disk space and memory for the `sequence` marker set (see `estimatelimits` in `config.yaml`)

## Quick start

```bash
# 1. Install the environment (Snakemake, Julia, bcftools, plink, ...)
pixi install

# 2. Set where the large raw tool output should go on your system
#    -> edit `tmpdir` in your profile, e.g. profiles/ghpcaudk/config.yaml

# 3. Check what would be run
pixi run dryrun

# 4. Small local test (2 replicates, array markers only)
pixi run test-local

# 5. Full run on the cluster (in the background, log in logs/snakemake/)
bash src/snakemake_ghpcaudk.sh                  # profiles/ghpcaudk
bash src/snakemake_ghpcaudk.sh profiles/genomedk
```

All software (the LAI tools, ADMIXTURE, Julia packages and the LocalAncestry.jl
system image) is downloaded and built automatically by the pipeline at the pinned
versions in `config.yaml`. Setup jobs run as cluster jobs, so they are compiled on
compute nodes, not on the login node. To only set up Julia:

```bash
pixi run snakemake --profile profiles/ghpcaudk setup_julia
```

Julia package caches and the system image are compiled for several x86-64 CPU
generations (`JULIA_CPU_TARGET` in `pixi.toml`), so they work on nodes with
different architectures. Install pixi itself on the login node with `pixi install`
as usual; conda packages are not CPU-specific.

Only one Snakemake instance can run in the repository at a time.

## Configuration

The experiment is set in [config.yaml](config.yaml); system-specific settings are in the profile.

| Key | Meaning |
|---|---|
| `replicates`, `generations`, `markersets` | Scenario grid |
| `lathresholds` | LocalAncestry.jl thresholds. The **first** one is used in tables and plots |
| `tools` | Which tools to run, e.g. remove `sparsepainter` to skip it |
| `run_block_comparison` | Run the LocalAncestry.jl block-size × threshold experiment |
| `software` | Pinned versions/commits of all tools |
| `estimatelimits` | Runtime (min) and memory (MB) requested per tool and marker set |

Override settings without editing the file:

```bash
pixi run snakemake --profile profiles/local --configfile config_test.yaml
# --config replaces the profile's `config:` entries, so repeat tmpdir
pixi run snakemake --profile profiles/local --config tmpdir=scratch 'replicates=[1,2,3]' 'tools=["localancestryjl","flare"]'
```

Prefer `--configfile` for experiment overrides. It does not interfere with the profile's `tmpdir`.

### Profiles

| Profile | Use | `tmpdir` |
|---|---|---|
| `profiles/local` | Run on the current machine | `scratch` (inside the repo) |
| `profiles/ghpcaudk` | SLURM via `cluster-generic` (writes job logs to `logs/slurm/`, which must exist) | ironbank research data |
| `profiles/genomedk` | SLURM via the Snakemake SLURM executor | `scratch` (change to your project storage) |

Each profile sets `tmpdir`, the location of the large raw tool output, in its `config:` section.
The pipeline stops with an error if no `tmpdir` is set.

All profiles set `rerun-triggers: mtime`. Jobs are only rerun when their input files
are newer than their output, not when rule code changes. Delete the outputs, or use
`--forcerun <rule>`, to rerun deliberately.

## Pipeline overview

| Step | Rule file | Output |
|---|---|---|
| 1. Data preparation | `rules/1. dataprep.smk` | Simulated data from Zenodo, 1000 Genomes reference, genomic maps, VCF subsets in `data/` |
| 2. Software setup | `rules/2. softwaresetup.smk` | Tool binaries, Julia packages and `deps/localancestry.so` |
| 3. Local ancestry | `rules/3. localancestry.smk` | Raw tool output in `<tmpdir>/localancestry/<tool>/`, runtime/memory in `benchmarks/<tool>/` |
| 4. Accuracy | `rules/4. benchmark.smk` | Per-locus accuracy in `results/accuracy/<tool>/`, summaries in `results/summary/` |
| 5. Block experiment | `rules/5. iablocksVSfixedblocks.smk` | `results/blocks/` |
| 6. Plots | `rules/6. plots.smk` | `plots/*.png` / `*.pdf` |

File names encode the scenario: `<replicate>_<markerset>_<g1>_<g2>[_<threshold>]`.

Accuracy is the fraction of target haplotypes whose inferred ancestry matches the
simulated ancestry, computed per locus. A *sequence* run is evaluated both on all loci
and on array loci only ("Sequence -> Array"). FLARE (clustered) and SparsePainter
report ancestries as indices. These are matched to population labels by the
permutation that maximizes overlap with the truth.

Draw the rule graph with `bash src/snakemake_plotdag.sh` (writes `plots/dag.png`).

### Real data (`rules/7. realdata.smk`)

LocalAncestry.jl is also run genome-wide (chromosomes 1–22) on an admixed 1000 Genomes
Phase 3 population (default: ACB). Results are compared with the published 1000 Genomes
local ancestry calls (RFMix v1.5.4,
[20140818_ancestry_deconvolution](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140818_ancestry_deconvolution/)).
Settings are under `realdata` in `config.yaml`, and `run_realdata` turns the analysis on or off.

- **Data:** the Phase 3 release VCFs (GRCh37, the same build as the published calls), restricted to biallelic SNPs with MAF ≥ 0.005.
- **References:** YRI (AFR) and CEU (EUR), the populations behind the published calls.
- **`totcm`:** each chromosome's genetic length, taken from the cM columns of the published calls.
- **Main metric:** diploid concordance, i.e. agreement in the number of AFR copies at each locus. It does not depend on phasing. Haplotype concordance is also reported, but the published "PopPhased" calls were re-phased, so their haplotype order can differ from the VCF.
- **Excluded loci:** loci where the published call is UNK (posterior < 0.9), NAT, centromere or miscall. Their share is reported as `published_excluded`.

Outputs in `results/realdata/`: `summary.csv` (per chromosome and genome-wide), `individuals.csv`
(global AFR proportion per individual, published vs. LocalAncestry.jl) and `resources.csv`
(runtime and peak memory per chromosome and genome-wide).

## Main outputs

- `results/summary/accuracy.csv`, `accuracy_sds.csv`: mean accuracy (and SD across replicates) per scenario and tool
- `results/summary/missing_accuracy.csv`: scenario/replicate combinations with no result
- `results/summary/{time,RAM}_use_{mean,sds}.csv`: runtime (minutes) and peak memory (GB)
- `results/summary/resources_table_latex.csv`: LaTeX-ready resource table
- `plots/accuracy_ancgen.png`, `accuracy_software.png`, `sensitivities.png`, `accuracy_vs_rfmix.png`, `accuracy_blocksize.png`

## Troubleshooting

- **Jobs fail immediately on the cluster with no log**: create `logs/slurm/` (the launcher script does this).
- **`deps/localancestry.so` missing**: run `pixi run snakemake --profile <cluster profile> setup_julia`.
- **Running accuracy scripts outside Snakemake**: `src/julia/prelim_accuracies.jl <tmpdir>` processes
  every finished tool output that has no accuracy file yet.
