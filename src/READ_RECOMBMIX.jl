
function read_recombmix(infile, outfile, chromosome)
    dump = readlines(infile)
    populations = replace.(split(replace(dump[1], r".+\: " => ""), '\t'), r"=.+" => "")
    popdict = Dict(string.(0:(length(populations)-1)) .=> populations)
    hapdict = Dict('0':'1' .=> 1:2)

    # Write
    file = open(outfile, "w")

    ## Write header
    println(file, "individual,chromosome,haplotype,basepairs,ancestry")

    # Write
    for s in dump[2:end]
        stmp = split(s, '\t')
        individual = replace(stmp[1], r"_.$" => "")
            for b in 4:3:length(stmp)
                println(file, join([individual, chromosome, hapdict[last(stmp[1])], stmp[b-2] * ":" * stmp[b-1], popdict[stmp[b]]], ","))

            end
    end
    close(file)
    return nothing
end

