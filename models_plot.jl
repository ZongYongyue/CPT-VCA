using CPTVCA
using ExactDiagonalization
using QuantumLattices
using Plots
using TimerOutputs
const vcatimer = TimerOutput()
#=
vca = loadData("line(1,2)U4_vca.jls")
#give a path in the reciprocal space 
k_path = ReciprocalPath(reciprocals(vca.unitcell.vectors), line"Γ-X", length=300)
#give the energy range
ω_range = range(-6, 6, length=400)
#caculate the spectral function of the VCA green function with half filling
G = singleParticleGreenFunction(vca, k_path, ω_range,Parameters(vca.origigenerator)[:U]/2;timer=vcatimer)
A = spectrum(G)
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function(L=2, U=4, n=1/2)",clims=(0, 3))
#savefig("tline(1,2)U4.pdf")
=#
#=
#square lattice t=-1, U=4 L=(2,2), n=1/2, NN
vca = loadData("squareL4U4_vca.jls")
k_path = ReciprocalPath(reciprocals(vca.unitcell.vectors), rectangle"Γ-X-M-Γ", length=100)
ω_range = range(-6, 6, length=400)
G = singleParticleGreenFunction(vca, k_path::ReciprocalPath,ω_range,Parameters(vca.origigenerator)[:U]/2;timer=vcatimer)
A = spectrum(G)
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function(L=2, U=4, n=1/2)",clims=(0, 3))
=#

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

fs = range(-3, 6, length=50)
rz = ReciprocalZone(reciprocals(vca.unitcell.vectors); length=300)
ef = energysurface(vca, rz, fs, Parameters(vca.origigenerator)[:U]/2, 2.7)
b = plot(rz, ef; color=:haline)
contour!(rz, m, levels=[0], linestyles=[:dash], linewidths=2, color=:white)

=#

#=
#square lattice, L = (3,4), n = 0.583
vca = loadData("squareL12_vca.jls")
path = @rectangle_str("Γ-X-M-Γ")
k_path = selectpath(BrillouinZone(reciprocals(vca.unitcell.vectors), (200, 200)), path)[1]
ω_range = range(-8, 8, length=400)
A = singleparticlespectrum(vca, k_path, ω_range, Parameters(vca.origigenerator)[:U]*(0.583); timer=vcatimer)
show(vcatimer)
#=
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="t1=-1, t2=0.3, t3=-0.2, U=4, L=(3,4), n=0.583",clims=(0, 3))
savefig("squareL12.pdf")
=#
=#