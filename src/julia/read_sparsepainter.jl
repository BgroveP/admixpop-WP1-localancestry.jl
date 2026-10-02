using DataFrames
using CodecZlib
using Base.Threads: @threads, ReentrantLock, @lock

function read_lines_any(path::AbstractString)
    io = endswith(path, ".gz") ? GzipDecompressorStream(open(path)) : open(path)
    lines = try
        readlines(io)
    finally
        close(io)
    end
    return lines
end

# argmax over a comma-separated probability string, e.g. "0.2,0.3,0.5" -> 3
@inline function anc_index(s::AbstractString)
    best = -1.0; idx = 0; k = 1
    for tok in eachsplit(s, ',')
        v = parse(Float64, tok)
        v > best && (best = v; idx = k)
        k += 1
    end
    return idx
end

const SegRow = NamedTuple{(:individual, :chromosome, :haplotype, :basepairs, :ancestry),
                          Tuple{String, String, Int, UnitRange{Int}, String}}

function read_sparsepainter(infile::AbstractString, chrom::AbstractString,
                            labels::Vector{<:AbstractString})
    labels = String.(labels)
    println("Reading lines from file")
    lines = read_lines_any(infile)

    data_start = findfirst(l -> !startswith(l, '#'), lines)
    data_start === nothing && error("No data lines found in $infile")

    header_parts = split(strip(lines[data_start]))
    header_parts[1] == "ind_name" ||
        error("Expected header starting with 'ind_name', got: $(lines[data_start])")
    positions = [parse(Int, p) for p in header_parts[2:end]]
    nbins = length(positions)

    rows = SegRow[]                       # concretely typed accumulator
    sizehint!(rows, 100_000)
    plock = ReentrantLock()

    println("Parsing lines to ancestry blocks")
    @threads for line in lines[(data_start + 1):end]
        (isempty(line) || all(isspace, line)) && continue
        parts = split(strip(line))
        ind = String(parts[1])
        cells = parts[2:end]
        nbins == length(cells) ||
            error("Row $ind has $(length(cells)) bins, expected $nbins")

        # Decode both haplotypes per bin into Int ancestry codes
        c1 = Vector{Int}(undef, nbins)
        c2 = Vector{Int}(undef, nbins)
        @inbounds for (i, cell) in enumerate(cells)
            bar = findfirst('|', cell)
            c1[i] = anc_index(SubString(cell, 1, bar - 1))
            c2[i] = anc_index(SubString(cell, bar + 1, length(cell)))
        end

        # Merge consecutive bins with the same ancestry into segments
        for (h, codes) in ((1, c1), (2, c2))
            prev = codes[1]; s0 = 1
            @inbounds for i in 2:nbins
                if codes[i] != prev            # Int comparison, no allocation
                    @lock plock push!(rows, (individual = ind, chromosome = chrom,
                                 haplotype = h,
                                 basepairs = positions[s0]:positions[i-1],
                                 ancestry = labels[prev]))
                    prev = codes[i]; s0 = i
                end
            end
            @lock plock push!(rows, (individual = ind, chromosome = chrom,
                         haplotype = h,
                         basepairs = positions[s0]:positions[nbins],
                         ancestry = labels[prev]))
        end
    end

    # Sort the output
    out = DataFrame(rows)
    sort!(out, [:individual, :chromosome, :haplotype, :basepairs])

    # Output
    return out   # already concretely typed, no column fixups needed
end