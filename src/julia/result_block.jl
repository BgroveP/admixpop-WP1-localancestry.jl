
using CSV
using DataFrames
using StatsBase
using AlgebraOfGraphics
using CairoMakie

# ── Load ─────────────────────────────────────────────────────────────
datadir = "results/blocks"
o = DataFrame[]
for f in readdir(datadir; join = true)
    thiso = CSV.read(f, DataFrame)
    parts = split(replace(replace(f, "." .=> "_"), "/" .=> "_"), "_")
    thiso.markerset .= parts[3]
    thiso.replicate .= parts[4]
    thiso.blocksize .= parse(Float64, parts[5] * "." * parts[6])
    thiso.IA       .= parse(Float64, parts[7] * "." * parts[8])
    push!(o, thiso)
end
o = vcat(o...)
deleteat!(o, .~o.array)
deleteat!(o, o.blocksize .== 0.0)
deleteat!(o, o.IA .< 0.4)

# ── Summarize ────────────────────────────────────────────────────────
acc_per_blocksize = combine(
    groupby(o, [:markerset, :blocksize, :IA]),
    :accuracy => mean  => :accuracy,
    :accuracy => length => :n, 
)
sort!(acc_per_blocksize, :blocksize)

# ── Translate ────────────────────────────────────────────────────────
acc_per_blocksize.markerset = replace.(string.(acc_per_blocksize.markerset), "array" => "Array", "sequence" => "Sequence")

# ── Plot ─────────────────────────────────────────────────────────────
set_theme!(theme_minimal(); fontsize = 16)   

plt = data(acc_per_blocksize) *
    mapping(
        :blocksize => "Max. block size",
        :accuracy  => "Accuracy",
        color = :IA => "IA criterion",
        group = :IA,
        col = :markerset => "Marker set",
    ) *
    (visual(Lines; linewidth = 2) +
     visual(Scatter; markersize = 6, strokecolor = :white, strokewidth = 0.5))

fig = draw(
    plt,
    figure = (; size = (1000, 450)),
    axis   = (; ygridcolor = :grey90, xgridvisible = false),
)

save("plots/accuracy_blocksize.png", fig; px_per_unit = 2)
save("plots/accuracy_blocksize.pdf", fig; px_per_unit = 2)


# 
combine(groupby(acc_per_blocksize, [:markerset, :blocksize]), :accuracy .=> [mean, minimum, median, maximum]) |> sort
sort(acc_per_blocksize, :accuracy)