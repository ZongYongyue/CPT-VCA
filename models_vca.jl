using CPTVCA
using ExactDiagonalization
using QuantumLattices
using LinearAlgebra
#=
#define a unitcell and a cluster
unitcell = Lattice([0, 0]; vectors=[[1, 0]])
cluster = Lattice([0, 0], [1, 0]; vectors = [[2, 0]])
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
#give a target space
target = TargetSpace(BinaryBases(4, 2), BinaryBases(4, 1), BinaryBases(4, 3))
#define the terms of Hamiltonian
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 4.0)
origiterms = (t, U)
t_r = Hopping(:t, -1.0, 1)
referterms = (t_r, U)
neighbors = Neighbors(0=>0.0, 1=>1.0)
#give the VCA method and save the vca data
vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, target; neighbors=neighbors, m=200)
saveData(vca, "line(1,2)U4_vca.jls")
=#

# square lattice, L = (2,2), n = 1/2
unitcell = Lattice([0, 0]; vectors=[[1, 0],[0, 1]])
cluster = Lattice(unitcell,(2,2),('p','p'))
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
target = TargetSpace(BinaryBases(8, 4), BinaryBases(8, 3), BinaryBases(8, 5))
t = Hopping(:t, Complex(-1.0), 1)
U = Hubbard(:U, Complex(4.0))
origiterms = (t, U)

function AFphase₁(bond::Bond)
    exp(im*dot([π,π],rcoordinate(bond)))
end
function AFphase₂(bond::Bond)
    -exp(im*dot([π,π],rcoordinate(bond)))
end
coupling₁ = Coupling(Index(:,FID{:f}(:,1//2,:)),Index(:,FID{:f}(:,1//2,:)))
coupling₂ = Coupling(Index(:,FID{:f}(:,-1//2,:)),Index(:,FID{:f}(:,-1//2,:)))
af₁ = Onsite(:af₁, Complex(1.0), coupling₁; amplitude = AFphase₁)
af₂ = Onsite(:af₂, Complex(1.0), coupling₂; amplitude = AFphase₂)

t_r = Hopping(:t, Complex(-1.0), 1)
referterms = (t_r, U, af₁, af₂)
neighbors = Neighbors(0=>0.0, 1=>1.0)
vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, target; neighbors=neighbors, m=200)
saveData(vca, "squareL4U4_vca.jls")

varparams = diag([(af₁ = a , af₂ = b) for a in range(-0.3,0.3,length=50), b in range(-0.3,0.3,length=50)])
vcas = VCAs(unitcell, cluster, hilbert, origiterms, referterms, target, varparams; neighbors=neighbors, m=200)
saveData(vcas, "squareL4U4_af.jls")
#=
# square lattice, L = (2,2), n = 3/4
unitcell = Lattice([0, 0]; vectors=[[1, 0],[0, 1]])
cluster = Lattice([0, 0], [1, 0],[0, 1], [1, 1]; vectors = [[2, 0], [0, 2]])
hilbert₁ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(unitcell))
hilbert₂ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
target = TargetSpace(BinaryBases(8, 6), BinaryBases(8, 5), BinaryBases(8, 7))
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 4.0)
origiterms = (t, U)
t_r = Hopping(:t, -1.0, 1)
referterms = (t_r, U)
neighbors = Neighbors(0=>0.0, 1=>1.0)
vca = VCA(unitcell, cluster, hilbert₁, hilbert₂, origiterms, referterms, target; neighbors=neighbors, m=100)
saveData(vca, "square(2,2)U4_vca.jls")
=#
#=
# square lattice, L = (2,3)
unitcell = Lattice([0, 0]; vectors=[[1, 0], [0, 1]])
cluster = Lattice([0, 0], [1, 0], [0, 1], [1, 1], [0, 2], [1, 2]; vectors = [[2, 0], [0, 3]])
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
target = TargetSpace(BinaryBases(12, 7), BinaryBases(12, 6), BinaryBases(12, 8))
t1 = Hopping(:t1, -1.0, 1)
t2 = Hopping(:t2, 0.4, 2)
U = Hubbard(:U, 8.0)
origiterms = (t1, t2, U)
t1_r = Hopping(:t1_r, -1.0, 1)
t2_r = Hopping(:t2_r, 0.4, 2)
referterms = (t1_r, t2_r, U)
neighbors = Neighbors(0=>0.0, 1=>1.0, 2=>√2)
@time vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, target; neighbors=neighbors, m=300)
#1.492265 seconds (5.24 M allocations: 327.429 MiB, 3.37% gc time, 97.45% compilation time: 38% of which was recompilation)
saveData(vca, "squareL6NNN_vca.jls")
=#
#=
# square lattice, L = (3,4), n = 0.583
unitcell = Lattice([0, 0]; vectors=[[1, 0], [0, 1]])
cluster = Lattice(unitcell, (3,4), ('p','p'))
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
target = TargetSpace(BinaryBases(24, 14), BinaryBases(24, 13), BinaryBases(24, 15))
t1 = Hopping(:t1, -1.0, 1)
t2 = Hopping(:t2, 0.3, 2)
t3 = Hopping(:t3, -0.2, 3)
U = Hubbard(:U, 4.0)
origiterms = (t1, t2, t3, U)
t1r = Hopping(:t1r, -1.0, 1)
t2r = Hopping(:t2r, 0.3, 2)
t3r = Hopping(:t3r, -0.2, 3)
referterms = (t1r, t2r,t3r, U)
neighbors = Neighbors(0=>0.0, 1=>1.0, 2=>√2, 3=>2.0)
@time vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, target; neighbors=neighbors, m=200)
#1008.520550 seconds (549.60 M allocations: 481.442 GiB, 1.93% gc time, 0.06% compilation time)
saveData(vca, "squareL12_vca.jls")
=#
#=
# honeycomb lattice L=6
unitcell = Lattice([0,0],[√3/2,1/2]; vectors = [[√3/2,3/2],[√3,0]])
cluster = Lattice([0,0],[√3/2,1/2],[√3,0],[√3,-1],[√3/2,-3/2],[0,-1]; vectors = [[3√3/2,3/2],[3√3/2,-3/2]])
hilbert₁ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(unitcell))
hilbert₂ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
k = EDKind(typeof((Hopping(:t, -1.0, 1),)))
table = Table(hilbert₂, Metric(EDKind{:FED}(), hilbert₂))
per = Perioder(unitcell, cluster, table)
=#
#=
#triangle lattice L=12
unitcell = Lattice([0, 0]; vectors = [[1, 0], [-1/2, √3/2]])
lattice = Lattice([0.0, 0.0], [1, 0], [-1/2, √3/2], [1/2, √3/2], [3/2, √3/2], 
[-1, √3], [0, √3], [1, √3], [2, √3], [-1/2, 3√3/2], [1/2, 3√3/2], [3/2, 3√3/2]; vectors = [[3, √3], [-3, √3]])
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
target = TargetSpace(BinaryBases(24, 12), BinaryBases(24, 11), BinaryBases(24, 13))
neighbors = Neighbors(0=>0.0, 1=>1.0, 2=>√3)
@time vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, target; neighbors=neighbors, m=200)
saveData(vca, "triangleL12U4_vca.jls")
=#