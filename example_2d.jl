using Base.Threads
using QuantumLattices
using ExactDiagonalization
using CPTVCA
using Plots
gr()
#launch multi-threads
Threads.nthreads() = 24
#define a 2d Lattice systerm as a cluster
lattice = Lattice([0, 0], [1, 0], [0, 1], [1, 1])
#give a bases
bases = BinaryBases(8, 4)
#define the terms of Hamiltonian
t= Hopping(:t, -1.0, 1)
U = Hubbard(:U, 4.0)
terms = (t, U)
neighbors = Neighbors(1=>1.0)
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))
#give the vectors to tile the space with the cluster 
supervectors = [[2, 0], [0, 2]]
#give the vectors of the unitcell
unitvectors = [[1, 0], [0, 1]]
#caculate the constant information needed in caculating the CPT green function 
info = CPTinfo(lattice, bases, terms, hilbert, neighbors, supervectors, unitvectors, (200, 200))
#give a path in the reciprocal space
k_path = @rectangle_str("Γ-X-M-Γ") 
#give the energy range
ω_range = range(-6, 6, length=400)
#caculate the spectral function of the CPT green function
A = CPTspec(info, k_path, ω_range)
#plot the the spectral function
specfig = heatmap(1:size(A)[2], ω_range, A, xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function of 2d Hubbard Model",clims=(0, 3))
N = size(A)[2]
ratios = [1, 1, sqrt(2)]
split_points = round.(Int, cumsum(ratios) ./ sum(ratios) * N)
special_points = pushfirst!(split_points,1)
tick_labels = ["Γ", "X", "M", "Γ"]
xticks!(special_points, tick_labels)
savefig("HubbardModel_2d.pdf")

