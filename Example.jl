using Arpack: eigs
using LinearAlgebra: Eigen
using SparseArrays: SparseMatrixCSC, spzeros
using QuantumLattices
using ExactDiagonalization

using Pkg;
Pkg.add("/Users/tzung/Library/Mobile Documents/com~apple~CloudDocs/mygit/CPTVCA")
using CPTVCA
#= 
"""
get  the matrix representation Hamiltonian of a single cluster systerm
"""
struct ClusterH{L<:AbstractLattice, G<:OperatorSum, M<:OperatorSum} <: Frontend
    cluster::L
    clH::G
    clHₘ::M
end
function ClusterH(cluster::AbstractLattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, targetspace::TargetSpace; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)
    ed = ED(cluster, hilbert, terms, targetspace; neighbors, boundary)
    clH = expand(ed.H)
    clHₘ = expand(ed.Hₘ)
    return ClusterH(cluster, clH, clHₘ)
end
=#
lattice = Lattice([0.0, 0.0], [1, 0])
bond = bonds(lattice, Neighbors(1=>1.0))
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))
bases = BinaryBases(4, 2)
coupling = Coupling(Index(:,FID{:f}(1,:,:)),Index(:,FID{:f}(1,:,:)))
t = Hopping(:t, Complex(-1.0), 1, coupling)
U = Hubbard(:U, Complex(4.0))
H = ClusterH(lattice, hilbert, (t, U), TargetSpace(bases); boundary = plain)

H.clH
H.clHₘ