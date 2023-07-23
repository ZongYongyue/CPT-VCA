using CPTVCA
using ExactDiagonalization
using QuantumLattices

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
vca = VCA(unitcell, cluster, hilbert, origiterms, referterms, target; neighbors=neighbors, m=100)
saveData(vca, "line(1,2)U4_vca.jls")

