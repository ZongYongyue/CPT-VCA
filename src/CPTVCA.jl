module CPTVCA

using Arpack: eigs
using LinearAlgebra: Eigen
using SparseArrays: SparseMatrixCSC, spzeros
using QuantumLattices
using ExactDiagonalization
using BlockArrays
export ClusterH, Cluster, Clusters
"""
define the cluster of a quantum lattice systerm
"""
struct Cluster{T<:Number} <: AbstractMatrix{T}
    intrasites::Matrix{T}
    id::Int
end
Base.size(c::Cluster)=size(c.intrasites)
Base.getindex(c::Cluster, i, j) = c.intrasites[i, j]
Base.setindex!(c::Cluster, v, i, j) = (c.intrasites[i, j] = v)
function Cluster(lattice::AbstractLattice, id::Int=0)
    intrasites = lattice.coordinates
    return Cluster(intrasites, id)
end

"""
get a supercluster consisting of a given cluster and its surrounding clusters 
"""
function Clusters(cluster::Cluster, steps::Tuple{Int, Vararg{Int}}, vectors::AbstractVector{<:AbstractVector{<:Number}})
    @assert length(steps) == length(vectors) "steps and vectors must have the same length"
    translations = Array{Array{Int64, 1}, 1}()
    if length(steps) == 2
        for i in -steps[1]:steps[1], j in -steps[2]:steps[2]
            push!(translations, [i, j])
        end
    elseif length(steps) == 3
        for i in -steps[1]:steps[1], j in -steps[2]:steps[2], k in -steps[3]:steps[3]
            push!(translations, [i, j, k])
        end
    else
        error("Dimension is limited to 2 or 3.")
    end
    supermatrix = tile(cluster, vectors, translations)
    index = convert(Int, size(supermatrix,2)/size(cluster.intrasites,2))
    supermatrix_block = BlockArray(supermatrix, [size(cluster.intrasites,1)],fill(size(cluster.intrasites,2),index))
    clusters = Vector{Cluster}()
    for i in 1:index
        push!(clusters, Cluster(supermatrix_block[Block(1,i)],i))
    end
    return clusters
end

"""
get the inter-cluster bonds of the given cluster and its surrounding clusters
"""
function interbonds(clusters::Vector{Cluster}, cluster::Cluster, neighbors::Neighbors)
    surroundings = filter(x->x!=cluster,clusters)
    for (index,value) in enumerate(surroundings)
        inter




"""
get the matrix representation Hamiltonian of a single cluster
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

end # module CPTVCA
