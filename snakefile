import re

# ============================================================
# Comparison of local ancestry inference tools
#   Default target: `finished.txt` (all summaries and plots).
#   See README.md for usage.
# ============================================================

configfile: "config.yaml"

# --- Variables from the config file (wildcard name in brackets) ---
replicates   = config["replicates"]     # {r}
generations  = config["generations"]    # {g}, {g1} (reference), {g2} (target)
populations  = config["populations"]    # {p}
cohorts      = config["cohorts"]        # {c}
markersets   = config["markersets"]     # {m}
lathresholds = config["lathresholds"]   # {t}
tools        = config["tools"]

ALL_MARKERSETS = ["array", "sequence"]  # genomic maps are always built for both
# Raw tool output location. Set per system in the profile (`config: ["tmpdir=..."]`)
# or on the command line with --config tmpdir=...
if not config.get("tmpdir"):
    raise ValueError(
        "No 'tmpdir' configured. Use a profile that sets it (e.g. --profile profiles/local) "
        "or add --config tmpdir=/path/to/large/storage. Note that passing --config on the "
        "command line replaces the profile's `config` entries, so include tmpdir=... as well."
    )
LA_DIR = str(config["tmpdir"]).rstrip("/") + "/localancestry"
DEFAULT_THRESHOLD = lathresholds[0]

KNOWN_TOOLS = ["localancestryjl", "flare", "flareclust", "rfmix", "recombmix", "sparsepainter", "gnomix"]
_unknown = set(tools) - set(KNOWN_TOOLS)
if _unknown:
    raise ValueError(f"Unknown tool(s) in config 'tools': {sorted(_unknown)}. Known tools: {KNOWN_TOOLS}")


def wc_constraint(values):
    alts = sorted((re.escape(str(v)) for v in values), key=len, reverse=True)
    return "|".join(alts)


wildcard_constraints:
    r = wc_constraint(replicates),
    g = wc_constraint(generations),
    g1 = wc_constraint(generations),
    g2 = wc_constraint(generations),
    p = wc_constraint(populations),
    c = wc_constraint(cohorts),
    m = wc_constraint(markersets),
    t = wc_constraint(lathresholds)


# --- Helpers ---
def _limits(tool, markerset):
    limits = config["estimatelimits"]
    return limits.get(tool, {}).get(markerset, limits["default"][markerset])


def _runtime(tool, markerset):
    return _limits(tool, markerset)["runtime"]


def _mem(tool, markerset):
    return _limits(tool, markerset)["mem_mb"]


def tool_files(pattern):
    """Expand `pattern` ({tool} plus the usual wildcards) for every configured tool.
    LocalAncestry.jl files additionally carry the threshold `_{t}`."""
    out = []
    for tool in tools:
        suffix = "_{t}" if tool == "localancestryjl" else ""
        out += expand(pattern.replace("{tool}", tool).replace("{suffix}", suffix),
                      r=replicates, m=markersets, g1=generations, g2=generations, t=lathresholds)
    return list(dict.fromkeys(out))  # de-duplicate, keep order


# Every Julia script needs the packages and system image created by `setup_julia`.
# `ancient` prevents a rebuild of the image from re-triggering the whole pipeline.
JULIA_ENV = ancient("deps/localancestry.so")


# --- Final target ---
def final_plots(wildcards):
    plots = ["plots/accuracy_ancgen.png", "plots/accuracy_software.png", "plots/sensitivities.png"]
    if config.get("run_block_comparison", True):
        plots.append("plots/accuracy_blocksize.png")
    if config.get("run_realdata", False):
        plots.append("results/realdata/summary.csv")
    return plots


localrules: comparison

rule comparison:
    input:
        final_plots
    output:
        "finished.txt"
    shell:
        "touch {output}"


include: "rules/1. dataprep.smk"
include: "rules/2. softwaresetup.smk"
include: "rules/3. localancestry.smk"
include: "rules/4. benchmark.smk"
include: "rules/5. iablocksVSfixedblocks.smk"
include: "rules/6. plots.smk"
include: "rules/7. realdata.smk"
