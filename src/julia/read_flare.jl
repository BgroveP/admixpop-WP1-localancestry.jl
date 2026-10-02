
using QGIO
using LocalAncestry

function read_flare(path)
    
    outcols = [:individual, :chromosome, :haplotype, :basepairs, :ancestry]

    populations = replace.(split(replace(QGIO.header(path)[end], r".+\<" => "", ">" => "", '\n' => ""), ","), r"=." => "")
    loci = QGIO.loci(path)
    samples = QGIO.samples(path)
    a = zeros(Int8, nrow(loci), 2 * length(samples))
    tmpancestries = zeros(Int8, length(samples))

    file = QGIO.open_vcf(path)
    buffer = QGIO.create_buffer()
    outdict = Dict(tuple.(loci.chromosome, loci.position, loci.identifier) .=> 1:nrow(loci))
    buffer[1:2] .= UInt8('#')
    #
    while buffer[1:2] == repeat([UInt8('#')], inner = 2)
        QGIO.readline!(file, buffer)
    end

    #
    while QGIO.readline!(file, buffer)
        outrow = outdict[(QGIO._buffer_chromosome(buffer), QGIO._buffer_position(buffer), QGIO._buffer_identifier(buffer))]
        _buffer_field!(tmpancestries,buffer, "AN1")
        a[outrow, 1:2:end] = tmpancestries .+ 1
        _buffer_field!(tmpancestries,buffer, "AN2")
        a[outrow, 2:2:end] = tmpancestries .+ 1
    end

    popdict = Dict(populations .=> 1)
    dicttranslate = Dict(keys(popdict) .=> populations)
    out = LocalAncestry.anc_reduce(a, samples, popdict, UnitRange.(loci.position, loci.position))
    out.basepairs = out.block
    out.chromosome .= loci.chromosome[1]

    for i in axes(out, 1)
        out.ancestry[i] = dicttranslate[out.ancestry[i]]
    end
    return out[:, outcols]
end


function _buffer_field!(c, b, f)
    # Initialize
    bufferposition = 1
    # Get necessary information
    fields = QGIO._buffer_field(b)
    focalfield = ifield(fields, f)
    # Scroll to data (one column before data )
    for _ in 1:8
        bufferposition = findnext(x -> x == UInt8('\t'), b, bufferposition) + 1
    end
    # Get haplotypes1
    if focalfield == 1
        for i in eachindex(c)
            bufferposition = findnext(x -> x == UInt8('\t'), b, bufferposition) + 1
            @views c[i] = b[bufferposition] - 48
            @views c[i+1] = b[bufferposition+2] - 48
        end
    elseif focalfield > 1
        for i in eachindex(c)
            bufferposition = findnext(x -> x == UInt8('\t'), b, bufferposition) + 1
            for _ in 1:(focalfield-1)
                bufferposition = findnext(x -> x == UInt8(':'), b, bufferposition) + 1
            end
            @views c[i] = b[bufferposition] - 48
        end
    else
        error("Did not find GT field in buffer")
    end
end

function nfields(f)

    return length(split(f, ":"))
end

function ifield(f, p)
    return findfirst(split(f, ":") .== p)
end

