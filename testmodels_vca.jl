using CPTVCA
using ExactDiagonalization
using QuantumLattices

# square lattice, L = (2,2)
unitcell = Lattice([0, 0]; vectors=[[1, 0],[0, 1]])
cluster = Lattice([0, 0], [1, 0],[0, 1], [1, 1]; vectors = [[2, 0], [0, 2]])
hilbert₁ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(unitcell))
hilbert₂ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
target = TargetSpace(BinaryBases(8, 4), BinaryBases(8, 3), BinaryBases(8, 5))
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 4.0)
origiterms = (t, U)
t_r = Hopping(:t, -1.0, 1)
referterms = (t_r, U)
neighbors = Neighbors(0=>0.0, 1=>1.0)
vca = VCA(unitcell, cluster, hilbert₁, hilbert₂, origiterms, referterms, target; neighbors=neighbors, m=100)
saveData(vca, "square(2,2)U4_vca.jls")




#= 
# honeycomb lattice L=6
unitcell = Lattice([0,0],[√3/2,1/2]; vectors = [[√3/2,3/2],[√3,0]])
cluster = Lattice([0,0],[√3/2,1/2],[√3,0],[√3,-1],[√3/2,-3/2],[0,-1]; vectors = [[3√3/2,3/2],[3√3/2,-3/2]])
hilbert₁ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(unitcell))
hilbert₂ = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
=#