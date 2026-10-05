
function read_vcf_field(path::AbstractString, c::Union{AbstractString,Int}, field::AbstractString)

    # Initialize
    loci = Vector{String}()
    individuals = Vector{String}()
    chromosome = ["", "chr"] .* string.([c, c])

    # Read n loci
    file = GZip.open(path)
    while !eof(file)
        try
            if any(readuntil(file, '\t') .== chromosome)
                push!(loci, readuntil(file, '\t'))
            end
            _ = readline(file)
        catch
            if !eof(file)
                error("Encountered unknown error while reading file")
            end
        end
    end

    # Read n individuals
    file = GZip.open(path)
    for line in eachline(file)
        if line[1:7] == "#CHROM\t"
            individuals = split(line, '\t')[10:end]
            break
        end
    end
    # Read numbers of individuals and loci
    nloci = length(loci)
    nindividuals = length(individuals)

    # Read haplotypes
    x = zeros(Int8, 2 * nindividuals, nloci)

    # Setup multi
    read_lock = ReentrantLock()
    write_lock = ReentrantLock()

    # For loop
    locidict = Dict{String,Int}(loci .=> 1:nloci)

    ## Scroll
    file = GZip.open(path)
    while readline(file)[1:7] != "#CHROM\t"
        Nothing
    end

    Threads.@threads for line in 1:nloci
        thisline = ""
        @lock read_lock thisline = readline(file)[1:(end-1)]
        split_string = split(thisline, '\t')

        if field == "GT"
            formatn = findfirst(split(split_string[9], ':') .== field)
            genotypes = getindex.(split.(split_string[10:end], ':'), formatn)
            @lock write_lock x[1:2:(2*nindividuals), locidict[split_string[2]]] = parse.(Int8, getindex.(genotypes, 1))
            @lock write_lock x[2:2:(2*nindividuals), locidict[split_string[2]]] = parse.(Int8, getindex.(genotypes, 3))
        elseif field == "AN"
            formatn = findfirst(split(split_string[9], ':') .== "AN1")
            genotypes = getindex.(split.(split_string[10:end], ':'), formatn)
            @lock write_lock x[1:2:(2*nindividuals), locidict[split_string[2]]] = parse.(Int8, getindex.(genotypes, 1))
            formatn = findfirst(split(split_string[9], ':') .== "AN2")
            genotypes = getindex.(split.(split_string[10:end], ':'), formatn)
            @lock write_lock x[2:2:(2*nindividuals), locidict[split_string[2]]] = parse.(Int8, getindex.(genotypes, 1))
        end
    end
    return x, individuals, loci
end

function read_vcf_header(path::AbstractString)

    header = Vector{String}()
    file = GZip.open(path)
    for line in eachline(file)
        if (line[1] == '#') & (line[2] == '#')
            push!(header, line)
        else
            break
        end
    end
    return header
end

function vcf_get_n_youngest(ga, n, group="target")

    if group == "target"
        use = ga.population .== "admixed"
    elseif (group == "ancestral") | (group == "current")
        use = (ga.population .!= "admixed") .& (length.(ga.population) .== length.(replace.(ga.population, r"anc$" => "")))
    else
        use = (length.(ga.population) .> length.(replace.(ga.population, r"anc$" => "")))

    end
    out = combine(groupby(ga[use, :], "population"), :individual => (x -> x[max(1, end - n + 1):end]) => :individual)
    return out
end
function slurm_isfinished(j)

end

function slurm_bulkinfo(user, date;
    info=["JobID", "Start", "JobName", "AllocCPUS", "State", "Elapsed", "MaxRSS"]
)
    s = split.(replace.(split(join(Char.(read(`sacct -u $(user) --format=$(info) --starttime=$(date) --units=K`))), '\n'), r" +" => " "), " ")

    out = DataFrame()
    for (i, c) in enumerate(info)
        out[:, c] = getindex.(s[3:(end-1)], i)
    end

    return out
end

function SLURM_wait(j)

    sleep(10)
    s = split(join(Char.(read(`sacct -j $(j)`))), "\n")

    # Wait until job starts
    while length(s) < 4
        sleep(4)
        s = split(join(Char.(read(`sacct -j $(j)`))), "\n")
    end

    # Wait until job finishes
    while all(split.(s, r" +")[3][6] .!= ["TIMEOUT", "COMPLETED"])
        sleep(4)
        s = split(join(Char.(read(`sacct -j $(j)`))), "\n")
    end
    return nothing
end

function SLURM_getmem(j)
    return parse(Float64, replace(replace(split(join(Char.(read(`sacct -j $(j) --format=MaxRSS --units=K`))), '\n')[4], " " => ""), "K" => ""))
end

function SLURM_gettime(j)

    hours, mins, secs = parse.(Int, split(split(join(Char.(read(`sacct -j $(j) --format=Elapsed`))), '\n')[4], ':'))
    return 60 * hours + mins + secs / 60
end

function SLURM_getstatus(j)

    return split.(replace.(split(join(Char.(read(`sacct -j $(j)`))), '\n'), r" +" => ' '), ' ')[4][5]
end


function SLURM_jobinfo(j, i)

    # Wait until the job is done
    while success(`squeue -j $(j)`)
        if (length(split(join(Char.(read(`squeue -j $(j)`))), '\n')) == 2)
            break
        end
        sleep(5)
    end

    # Get info
    if i == "memory"
        return SLURM_getmem(j)
    elseif i == "time"
        return SLURM_gettime(j)
    elseif i == "status"
        return SLURM_getstatus(j)
    end
end

function evaluate(predla, truela, gmap)

    nhaplotypes = nrow(unique(predla[:, ["individual", "haplotype"]]))
    outdict = Dict(gmap[!, "position(bp)"] .=> 1:nrow(gmap))
    out = DataFrame(locus=1:nrow(gmap), ncorrect=0, incorrect=0, unassigned=nhaplotypes)
    wlock = ReentrantLock()
    @threads for i in axes(truela, 1)
        tmp = predla[(predla.individual.==truela.individual[i]).&(predla.haplotype.==truela.haplotype[i]).&(length.(intersect.([truela.basepairs[i]], predla.basepairs)).>0), :]

        if nrow(tmp) > 0
            itstart = searchsortedfirst(gmap[!, "position(bp)"], first(truela.basepairs[i]))
            tstart = gmap[itstart, "position(bp)"]
            tend = gmap[max(itstart, searchsortedfirst(gmap[!, "position(bp)"], last(truela.basepairs[i])) - 1), "position(bp)"]
            for j in axes(tmp, 1)
                start = outdict[max(first(tmp.basepairs[j]), tstart)]
                stop = outdict[min(last(tmp.basepairs[j]), tend)]

                @lock wlock out.unassigned[start:stop] .-= 1

                if truela.ancestry[i] == tmp.ancestry[j]
                    @lock wlock out.ncorrect[start:stop] .+= 1
                else
                    @lock wlock out.incorrect[start:stop] .+= 1
                end
            end
        end
    end
    out.pcorrect = out.ncorrect ./ nhaplotypes
    out.pincorrect = out.incorrect ./ nhaplotypes
    out.punassigned = out.unassigned ./ nhaplotypes
    return out
end


function expand_grid(; kwargs...)
    thistuple = Base.product(values(kwargs)...) |> collect |> pairs
    out = DataFrame()
    for (i, c) in enumerate(keys(kwargs))
        out[:, c] = [thistuple[x][i] for x in 1:length(thistuple)]
    end
    sort!(out)
    return out
end


function longest(x)
    return x[length.(x).==maximum(length.(x))]
end
function longest_row(x)
    return length.(x).==maximum(length.(x))
end

function parse_time(x, unit="minutes")
    s = parse.(Int, split(replace(x, "-" => ":"), ":"))
    scales = [24 * 60, 60, 1, 1 / 60]
    if length(s) == 4
        o = sum(s .* scales)
    elseif length(s) == 3
        o = sum(s .* scales[2:end])
    else
        error("wrong format: " * x)
    end
    return o
end

function parse_memory(x)
    scale = 1.0
    if any(x[end] .== ['K', 'M', 'G'])
        if x[end] == 'K'
            scale = 1000
        elseif x[end] == 'M'
            scale = 10^6
        else
            scale = 10^9
        end
        return parse(Float64, x[1:(end-1)])
    else
        return parse(Float64, x)
    end
end

function overlap_matrix(tru, est)
    joined = innerjoin(tru, est, on = [:individual, :haplotype, :chromosome],
                       makeunique = true, validate = (false, false))
    labels1 = unique(tru.ancestry)
    labels2 = unique(est.ancestry)

    M = zeros(Int, length(labels1), length(labels2))
    for r in eachrow(joined)
        ov = length(intersect(r.basepairs, r.basepairs_1))
        ov == 0 && continue
        j = findfirst(==(r.ancestry),   labels1)   # true label
        i = findfirst(==(r.ancestry_1), labels2)   # est label
        M[i, j] += ov
    end
    return (M, labels1, labels2)
end
