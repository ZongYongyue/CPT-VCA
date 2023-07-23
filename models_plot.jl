using CPTVCA
using ExactDiagonalization
using QuantumLattices
using Plots
using TimerOutputs
const vcatimer = TimerOutput()

#square lattice t=-1, U=4 L=(2,2), n=1/2, NN
vca = loadData("square(2,3)U8_vca.jls")
path = @rectangle_str("Γ-X-M-Γ")
k_path = selectpath(BrillouinZone(reciprocals(vca.unitcell.vectors), (200, 200)), path)[1]
ω_range = range(-8, 8, length=400)
A = singleparticlespectrum(vca, k_path, ω_range, Parameters(vca.origigenerator)[:U]/2; timer=vcatimer)
show(vcatimer)
#=
──────────────────────────────────────────────────────────────────────
Time                    Allocations      
───────────────────────   ────────────────────────
Tot / % measured:       90.4s /  98.4%           94.6GiB /  99.6%    

Section     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────
spectrum         1    88.9s  100.0%   88.9s   94.3GiB  100.0%  94.3GiB
GFpathv        1    67.2s   75.6%   67.2s   81.5GiB   86.5%  81.5GiB
CGFv           1    21.7s   24.4%   21.7s   12.8GiB   13.5%  12.8GiB
──────────────────────────────────────────────────────────────────────
=#
#f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="t=-1, U=0, L=(2,2), n=1/2, NN",clims=(0, 3))
#savefig("square(2,3)U8.pdf")

#=
#square lattice t=-1, U=4 L=(2,2), n=3/4, NN
vca = loadData("square(2,2)U4_vca.jls")
path = (points = ((1//4, 1//4), (1//2, 0//1), (1//2, 1//2), (0//1, 0//1), (1//2, 0//1)), labels = ("H", "X", "M", "Γ", "X"))
k_path = selectpath(BrillouinZone(reciprocals(vca.unitcell.vectors), (200, 200)), path)[1]
ω_range = range(-6, 6, length=300)
@time A = singleparticlespectrum(vca, k_path, ω_range, Parameters(vca.origigenerator)[:U]*(3/4))
# 78.434298 seconds (405.20 M allocations: 31.885 GiB, 3.53% gc time, 0.21% compilation time)
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="t=-1, U=4, L=(2,2), n=3/4, NN",clims=(0, 3))
savefig("square(2,2)U4n3to4.pdf")
=#
#=
#square lattice t=-1 t'=0.4, U=8, L=(2,3), n=2/3, NN,NNN
vca = loadData("squareL6NNN_vca.jls")
#path = @rectangle_str("Γ-X-M-Γ")
path = (points = ((1//4, 1//4), (1//2, 0//1), (1//2, 1//2), (0//1, 0//1), (1//2, 0//1)), labels = ("H", "X", "M", "Γ", "X"))
k_path = selectpath(BrillouinZone(reciprocals(vca.unitcell.vectors), (200, 200)), path)[1]
ω_range = range(-8, 8, length=400)
@time A = singleparticlespectrum(vca, k_path, ω_range, Parameters(vca.origigenerator)[:U]*(2/3))
#337.968639 seconds (2.49 G allocations: 160.181 GiB, 4.31% gc time, 0.05% compilation time)
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="t1=-1, t2=0.4, U=8, L=(2,3), n=2/3, NN,NNN",clims=(0, 3))
savefig("tsquareL6NNN.pdf")
=#
#=
#square lattice, L = (3,4), n = 0.583
vca = loadData("squareL12NNNN_vca.jls")
path = @rectangle_str("Γ-X-M-Γ")
k_path = selectpath(BrillouinZone(reciprocals(vca.unitcell.vectors), (200, 200)), path)[1]
ω_range = range(-8, 8, length=400)
@time A = singleparticlespectrum(vca, k_path, ω_range, Parameters(vca.origigenerator)[:U]*(0.583))
#701.832146 seconds (8.64 G allocations: 475.984 GiB, 6.11% gc time, 0.04% compilation time)
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="t1=-1, t2=0.3, t3=-0.2, U=4, L=(3,4), n=0.583",clims=(0, 3))
savefig("squareL12NNNN.pdf")
=#