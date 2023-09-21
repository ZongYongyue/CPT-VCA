using CPTVCA
using ExactDiagonalization
using QuantumLattices
using KrylovKit
using Arpack
using LinearAlgebra
#=
#square lattice
unitcell = Lattice([0, 0]; vectors=[[1, 0],[0, 1]])
cluster = Lattice(unitcell,(2,2),('p','p'))
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
cons= Conservation(4, 1//2; np=4, sz=0)
t = Hopping(:t, Complex(-1.0), 1)
U = Hubbard(:U, Complex(8.0))
μ = Onsite(:μ, Complex(-4.0))
origiterms = (t, U, μ)
t_r = Hopping(:t, Complex(-1.0), 1)
#af = Onsite(:af, Complex(0.0), MatrixCoupling(:, FID, :, σ"z", :); amplitude=antiferromagnetism([π, π]))
sp1 = Onsite(:sp1, Complex(0.0), MatrixCoupling(:, FID, :, σ"z", :); amplitude=belongs([1,4]))
sp2 = Onsite(:sp2, Complex(0.0), (-1)*MatrixCoupling(:, FID, :, σ"z", :); amplitude=belongs([2,3]))
#referterms = (t_r, U, μ, af)
referterms = (t_r, U, μ, sp1, sp2)
neighbors = Neighbors(0=>0.0, 1=>1.0)
=#

#triangle lattice
unitcell = Lattice([0, 0]; vectors=[[1, 0],[-1/2, √3/2]])
cluster = Lattice([0, 0],[1, 0], [1/2, √3/2];vectors=[[3/2,√3/2],[0,√3]])
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
cons = Conservation(3,1//2;np=3)
t = Hopping(:t, Complex(-1.0), 1)
U = Hubbard(:U, Complex(12.0))
μ = Onsite(:μ, Complex(-6.0))
origiterms = (t, U, μ)

r = Onsite(:r, Complex(0.0),spinrotation(π/2, 7π/6); amplitude=belongs([1]))
g = Onsite(:g, Complex(0.0),spinrotation(π/2, 11π/6); amplitude=belongs([2]))
b = Onsite(:b, Complex(0.0),spinrotation(π/2, π/2); amplitude=belongs([3]))

referterms = (t, U, r, g, b, μ)
neighbors = Neighbors(0=>0.0, 1=>1.0)



varparams = [(r = a, g = a, b = a) for a in range(0,0.3,50)]
#varparams = [(sp1 = a, sp2 = a) for a in range(0,0.3,50)]
rz = ReciprocalZone(reciprocals(cluster.vectors); length=100)

#=
vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, cons; neighbors=neighbors, m=200)
k_path = ReciprocalPath(reciprocals(vca.unitcell.vectors), hexagon"Γ-M₂-K-Γ, 120°", length=100)
ω_range = range(-8, 8, length=400)
fq = ω_range .+ (0.05*im)
G = singleParticleGreenFunction(vca, k_path, fq)
A = spectrum(G)
f = plot(k_path, ω_range, A; xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function",clims=(0, 3))
=#


