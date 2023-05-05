module ClusterPerturbationT

using Arpack: eigs
using LinearAlgebra: Eigen
using SparseArrays: SparseMatrixCSC, spzeros
using QuantumLattices
using ExactDiagonalization

export ClusterH

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










end #module