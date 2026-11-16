

# Standalone accuracy runner — bypasses Snakemake (no locks, no metadata).
# Sequential: one job at a time, skips existing outputs and unfinished inputs.

using Dates
using Base.Threads

const TMPDIR   = "/usr/home/ironbank/researchdata/bgp_public/admixpop-WP1-localancestry.jl"
const DATA_DIR = "data"
const OUT_ROOT = "results/accuracy"
const MIN_AGE  = Minute(5)

# (est subdir, script, output subdir, glob pattern, name function, threads)
methods = [
    ("flare",           "src/julia/accuracy_flare.jl",           "flare",           "*.anc.vcf.gz",                      f -> replace(basename(f), r"\.anc\.vcf\.gz$" => ""), 1),
    ("flareclust",      "src/julia/accuracy_flareclust.jl",      "flareclust",      "*.anc.vcf.gz",                      f -> replace(basename(f), r"\.anc\.vcf\.gz$" => ""), 1),
    ("rfmix",           "src/julia/accuracy_rfmix.jl",           "rfmix",           "*.msp.tsv",                         f -> replace(basename(f), r"\.msp\.tsv$" => ""), 1),
    ("recombmix",       "src/julia/accuracy_recombmix.jl",       "recombmix",       "*/inferred_local_ancestral_values.txt", f -> basename(dirname(f)), 1),
    ("gnomix",          "src/julia/accuracy_gnomix.jl",          "gnomix",          "*/query_results.msp",               f -> basename(dirname(f)), 1),
    ("localancestryjl", "src/julia/accuracy_localancestryjl.jl", "localancestryjl", "*.csv",                             f -> replace(basename(f), r"\.csv$" => ""), 1),
    ("sparsepainter",   "src/julia/accuracy_sparsepainter.jl",   "sparsepainter",   "*_prob.txt.gz",                     f -> replace(basename(f), r"_prob\.txt\.gz$" => ""), 10),
]

using Dates

is_stable(p) = isfile(p) && (now(UTC) - unix2datetime(mtime(p))) > MIN_AGE
truth_file(name) = joinpath(DATA_DIR, "replicate_$(first(split(name, '_')))", "local_ancestries.csv")

function find_jobs()
    jobs = NamedTuple[]
    skipped = 0
    for (sub, script, outdir, pattern, namefn, nt) in methods
        base = joinpath(TMPDIR, "localancestry", sub)
        isdir(base) || continue
        suffix = startswith(pattern, "*/") ? split(pattern, '/', limit = 2)[2] : pattern[2:end]
        for entry in readdir(base)
            est = joinpath(base, entry)
            if startswith(pattern, "*/")
                isdir(est) || continue
                inner = joinpath(est, suffix)
                isfile(inner) || continue
                est = inner
            else
                isfile(est) && endswith(entry, suffix) || continue
            end
            name = namefn(est)
            out  = joinpath(OUT_ROOT, outdir, "$name.csv")
            tru  = truth_file(name)
            if isfile(out) || !is_stable(est) || !isfile(tru)
                skipped += 1
                continue
            end
            push!(jobs, (est = est, tru = tru, out = out, script = script, nt = nt))
        end
    end
    sort!(jobs, by = j -> j.est)
    jobs, skipped
end

find_jobs()[1]  |> length

function main()
    jobs, skipped = find_jobs()
    println("Found $(length(jobs)) jobs to run, $skipped skipped (existing/unfinished).\n")
    failed = String[]
    @threads for i in rand(axes(jobs, 1), length(jobs))
        j = jobs[i]
        mkpath(dirname(j.out))
        println("[$i/$(length(jobs))] $(j.script) $(basename(j.est))")
        cmd = j.nt > 1 ? `julia -t $(j.nt) --project=. $(j.script) $(j.est) $(j.tru) $(j.out)` : `julia --project=. $(j.script) $(j.est) $(j.tru) $(j.out)`

        if !isfile(j.out)
            try
                success(run(cmd)) || push!(failed, j.est)
                println("DONE [$i/$(length(jobs))] $(j.script) $(basename(j.est))")
            catch e
                println(stderr, "✗ failed: $(j.est) — $(typeof(e))")
                push!(failed, j.est)
            end
        end
    end
    println("\nDone. $(length(jobs) - length(failed)) succeeded, $(length(failed)) failed.")
    isempty(failed) || foreach(f -> println("  ✗ $f"), failed)
end

main()
