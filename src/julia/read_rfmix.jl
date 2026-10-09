
function read_rfmix(file::String, loci)

    # Constants
    outcols = [:individual, :chromosome, :haplotype, :basepairs, :ancestry]

    # Reads
    ancestry = CSV.read(file, DataFrame, header=2, types=String)
    ancestry = stack(ancestry, 7:ncol(ancestry))
    populations = replace.(split(replace(readline(file), r".+: " => ""), "\t"), r"=.+" .=> [""])
    enddict = Dict{String,Int}(string.(loci.position[2:end]) .=> loci.position[1:(end-1)])

    ancestry[:, "chromosome"] = ancestry[:, "#chm"]
    ancestry[:, "individual"] = replace.(ancestry.variable, r"\..+" => "")
    ancestry[:, "haplotype"] = parse.(Int, replace.(ancestry.variable, r".+\." => "")) .+ 1
    ancestry[:, "basepairs"] = UnitRange.(parse.(Int, ancestry.spos), get.(Ref(enddict), ancestry.epos, missing))
    ancestry[:, "ancestry"] = get.(Ref(populations), parse.(Int, ancestry.value) .+ 1, missing)

    # Remove consecutive samples
    ancestry.startrow = axes(ancestry, 1)
    ancestry.endrow = axes(ancestry, 1)
    ancestry.discard .= false

    for i in axes(ancestry, 1)[2:end]
        if ancestry[i-1, outcols[[1, 2, 3, 5]]] == ancestry[i, outcols[[1, 2, 3, 5]]]

            ancestry.discard[i] = true
            ancestry.startrow[i] = ancestry.startrow[i-1]
            ancestry.endrow[ancestry.startrow[i]] = ancestry.endrow[i]
        end
    end


    ancestry.basepairs = UnitRange.(first.(ancestry.basepairs[ancestry.startrow]), last.(ancestry.basepairs[ancestry.endrow]))
    deleteat!(ancestry, ancestry.discard)

    #
    return ancestry[:, outcols]
end

function read_rfmix2(file::String, loci)
    pops = map(p -> replace(p, r"=.*" => ""),
               split(replace(readline(file), r".+: " => ""), '\t'))

    pos = loci.position
    enddict = Dict{eltype(pos),eltype(pos)}(pos[2:end] .=> pos[1:end-1])

    df = CSV.read(file, DataFrame; header=2)
    chm  = df[!, "#chm"]
    spos = df[!, "spos"]
    epos = df[!, "epos"]
    samplecols = names(df)[7:end]
    n = nrow(df)

    ind = String[];  sizehint!(ind, length(samplecols) * n ÷ 4)
    chr = Vector{eltype(chm)}();  hap = Int[]
    bps = Vector{UnitRange{eltype(pos)}}()
    anc = Vector{Union{String,Missing}}()

    for name in samplecols
        indiv = replace(name, r"\..+" => "")
        hapn  = parse(Int, replace(name, r".+\." => "")) + 1
        vals  = df[!, name]
        i = 1
        while i <= n
            v = vals[i]
            j = i
            while j < n && vals[j+1] == v   # extend the run in place
                j += 1
            end
            push!(ind, indiv)
            push!(chr, chm[i])
            push!(hap, hapn)
            push!(bps, UnitRange(spos[i], get(enddict, epos[j], missing)))
            push!(anc, get(pops, v + 1, missing))
            i = j + 1
        end
    end

    return DataFrame(individual = ind, chromosome = chr,
                     haplotype = hap, basepairs = bps, ancestry = anc)
end