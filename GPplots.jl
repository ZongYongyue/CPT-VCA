using Plots
using CPTVCA

M = 0.3
L = 101
color = palette(:Spectral_4)
#=
gps2 = loadData("./temp/squareL4U2af_gp.jls")
gps4 = loadData("./temp/squareL4U4af_gp.jls")
gps8 = loadData("./temp/squareL4U8af_gp.jls")
gps16 = loadData("./temp/squareL4U16af_gp.jls")

max_val = maximum(gps2)
Δgps2 = gps2 .- max_val
yticks = [-0.10, -0.08, -0.06, -0.04, -0.02, 0.0]
ymin, ymax = -0.101, 0.001
plot(range(-M, M, L), Δgps2, linecolor=color[1], label=nothing, yticks=yticks, ylim=(ymin, ymax), xlabel="M", ylabel="Ω-Ω₀", title=" Néel Weiss field , half-filling, 2×2 cluster")
min_val = minimum(Δgps2)
min_idx = argmin(Δgps2)
scatter!(range(-M, M, L), Δgps2, markercolor=color[1], markersize=4, label="U=2")
scatter!([range(-M, M, L)[101-min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps4)
Δgps4 = gps4 .- max_val
plot!(range(-M, M, L), Δgps4, linecolor=color[2], label=nothing)
min_val = minimum(Δgps4)
min_idx = argmin(Δgps4)
scatter!(range(-M, M, L), Δgps4, markercolor=color[2], markersize=4, label="U=4")
scatter!([range(-M, M, L)[101-min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps8)
Δgps8 = gps8 .- max_val
plot!(range(-M, M, L), Δgps8, linecolor=color[3], label=nothing)
min_val = minimum(Δgps8)
min_idx = argmin(Δgps8)
scatter!(range(-M, M, L), Δgps8, markercolor=color[3], markersize=4, label="U=8")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps16)
Δgps16 = gps16 .- max_val
plot!(range(-M, M, L), Δgps16, linecolor=color[4], label=nothing)
min_val = minimum(Δgps16)
min_idx = argmin(Δgps16)
scatter!(range(-M, M, L), Δgps16, markercolor=color[4], markersize=4, label="U=16")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label="minimum")
=#

#=
gps4 = loadData("./temp/triangleL3U4af_gp.jls")
gps5 = loadData("./temp/triangleL3U5af_gp.jls")
gps8 = loadData("./temp/triangleL3U8af_gp.jls")
gps12 = loadData("./temp/triangleL3U12af_gp.jls")


max_val = maximum(gps4)
Δgps4 = gps4 .- max_val
plot(range(-M, M, L), Δgps4, linecolor=color[1], label=nothing, xlabel="M", ylabel="Ω-Ω₀", title=" 120⁰ Néel, half-filling, triangle Lattice(L=3),ϕ=0")
min_val = minimum(Δgps4)
min_idx = argmin(Δgps4)
scatter!(range(-M, M, L), Δgps4, markercolor=color[1], markersize=4, label="U=4")
scatter!([range(-M, M, L)[101-min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps5)
Δgps5 = gps5 .- max_val
plot!(range(-M, M, L), Δgps5, linecolor=color[2], label=nothing)
min_val = minimum(Δgps5)
min_idx = argmin(Δgps5)
scatter!(range(-M, M, L), Δgps5, markercolor=color[2], markersize=4, label="U=5")
scatter!([range(-M, M, L)[101-min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps8)
Δgps8 = gps8 .- max_val
plot!(range(-M, M, L), Δgps8, linecolor=color[3], label=nothing)
min_val = minimum(Δgps8)
min_idx = argmin(Δgps8)
scatter!(range(-M, M, L), Δgps8, markercolor=color[3], markersize=4, label="U=8")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps12)
Δgps12 = gps12 .- max_val
yticks = [-0.10, -0.08, -0.06, -0.04, -0.02, 0.0]
ymin, ymax = -0.101, 0.001
plot!(range(-M, M, L), Δgps12, linecolor=color[4], label=nothing)
min_val = minimum(Δgps12)
min_idx = argmin(Δgps12)
scatter!(range(-M, M, L), Δgps12, markercolor=color[4], markersize=4, label="U=12")
scatter!([range(-M, M, L)[101-min_idx]], [min_val], markercolor=:black, markersize=7, label="minimum")
=#
#=
gps4 = loadData("./temp/triangleL3U4Mπ6thaf_gp.jls")
gps5 = loadData("./temp/triangleL3U5Mπ6thaf_gp.jls")
gps8 = loadData("./temp/triangleL3U8Mπ6thaf_gp.jls")
gps12 = loadData("./temp/triangleL3U12Mπ6thaf_gp.jls")
min = 51
L = 100
splice!(gps4, min)
splice!(gps5, min)
splice!(gps8, min)
splice!(gps12, min)
max_val = maximum(gps4)
Δgps4 = gps4 .- max_val
plot(range(-M, M, L), Δgps4, linecolor=color[1], label=nothing, xlabel="M", ylabel="Ω-Ω₀", title=" 120⁰ Néel, half-filling, triangle Lattice(L=3),ϕ=π/6")
min_val = minimum(Δgps4)
min_idx = argmin(Δgps4)
scatter!(range(-M, M, L), Δgps4, markercolor=color[1], markersize=4, label="U=4")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps5)
Δgps5 = gps5 .- max_val
plot!(range(-M, M, L), Δgps5, linecolor=color[2], label=nothing)
min_val = minimum(Δgps5)
min_idx = argmin(Δgps5)
scatter!(range(-M, M, L), Δgps5, markercolor=color[2], markersize=4, label="U=5")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps8)
Δgps8 = gps8 .- max_val
plot!(range(-M, M, L), Δgps8, linecolor=color[3], label=nothing)
min_val = minimum(Δgps8)
min_idx = argmin(Δgps8)
scatter!(range(-M, M, L), Δgps8, markercolor=color[3], markersize=4, label="U=8")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label=nothing)

max_val = maximum(gps12)
Δgps12 = gps12 .- max_val
plot!(range(-M, M, L), Δgps12, linecolor=color[4], label=nothing)
min_val = minimum(Δgps12)
min_idx = argmin(Δgps12)
scatter!(range(-M, M, L), Δgps12, markercolor=color[4], markersize=4, label="U=12")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label="minimum")
=#

max_val = maximum(gps)
Δgps = gps .- max_val
plot(range(0, M, L), Δgps, linecolor=color[4], label=nothing)
min_val = minimum(Δgps)
min_idx = argmin(Δgps)
scatter!(range(-M, M, L), Δgps, markercolor=color[4], markersize=4, label="U=8")
scatter!([range(-M, M, L)[min_idx]], [min_val], markercolor=:black, markersize=7, label="minimum")
