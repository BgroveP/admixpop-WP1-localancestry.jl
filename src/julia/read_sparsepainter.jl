
using DataFrames
using CodecZlib

# Read a plain or .gz-compressed text file into a Vector{String}
function read_lines_any(path::AbstractString)
    io = endswith(path, ".gz") ? GzipDecompressorStream(open(path)) : open(path)
    lines = try
        readlines(io)
    finally
        close(io)
    end
    return lines
end

function parse_onehot(cell::Vector{Float64}, labels::Vector{String})
    # "1,0,0" -> index of the nonzero entry (one-hot), or error
    _, idx = findmax(cell)
    return labels[idx]
end

function read_sparsepainter(infile::AbstractString, chrom::AbstractString, labels::Vector{<:AbstractString})
    
    # Read lines from file
    lines = read_lines_any(infile)

    # Skip any leading comment lines (e.g. "#Storage Mode: Raw")
    data_start = findfirst(l -> !startswith(strip(l), '#'), lines)
    data_start === nothing && error("No data lines found in $infile")

    # Header line: "ind_name <pos> <pos> ..."
    header_parts = split(strip(lines[data_start]))
    header_parts[1] == "ind_name" || error("Expected header starting with 'ind_name', got: $(lines[data_start])")
    positions = parse.(Int, header_parts[2:end])
    
    rows = NamedTuple[]
    for line in lines[(data_start + 1):end]
        strip(line) == "" && continue
        parts = split(strip(line))
        ind = parts[1]
        cells = parts[2:end]
        nbins = length(cells)
        haps = parse.(Float64, reduce(hcat, split.(reduce(hcat, split.(cells, '|')), ',')))
        iranges = [1:2:(2*length(positions)), 2:2:(2*length(positions))]

        for h in (1, 2)
            prev_anc  = ""
            seg_start = 0
            for i in 1:nbins
                anc = parse_onehot(haps[:, iranges[h][i]], labels)
                if prev_anc == ""
                    seg_start = positions[i]
                    prev_anc = anc
                elseif prev_anc != anc
                    push!(rows, (individual = ind, 
                                chromosome = chrom,
                                haplotype = h, 
                                basepairs = seg_start:positions[i-1],
                                ancestry = prev_anc))
                    prev_anc = anc
                    seg_start = positions[i]
                end
            end
            push!(rows, (individual = ind, 
                        chromosome = chrom,
                        haplotype = h, 
                        basepairs = seg_start:positions[i-1],
                        ancestry = prev_anc)
                        )
                    end
    end

    est = DataFrame(rows)
    # Column types to match the target format exactly
    est.individual  = String.(est.individual)
    est.chromosome  = String.(est.chromosome)
    est.haplotype   = Int.(est.haplotype)
    est.ancestry    = String.(est.ancestry)

    return est
end

est = main()