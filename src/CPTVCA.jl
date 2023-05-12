module CPTVCA

using Arpack: eigs
using LinearAlgebra: Eigen, dot
using SparseArrays: SparseMatrixCSC, spzeros
using QuantumLattices
using ExactDiagonalization
using BlockArrays
export Cluster, Clusters, interhopping, ClusterH
"""
define the cluster of a quantum lattice systerm
"""
struct Cluster{T<:Number} <: AbstractMatrix{T}
    intrasites::Matrix{T}
    id::Int
end
Base.size(c::Cluster)=size(c.intrasites)
Base.length(c::Cluster)=size(c.intrasites,2)
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
obtain the inter-cluster hopping  between the given cluster and its surrounding clusters with specific momentum k : V^{0,m}(k) 
"""
function interhopping(clusters::Vector{Cluster}, cluster::Cluster, hilbert::Hilbert, neighbors::Neighbors, k::Vector)
    surroundings = filter(x->x!=cluster,clusters)
    spin = Rational((hilbert[1].nspin-1)/2)
    norbital = hilbert[1].norbital
    ops = Vector{Operator}()
    for (index,value) in enumerate(surroundings), s in -spin:spin, n in 1:norbital
        seqs = interlinks(value, cluster, neighbors)
        if !isempty(seqs)
            for q in 1:length(seqs)
                op = exp(-im*dot(k, value[:,1]))*CompositeIndex(Index(seqs[q][2],FID{:f}(n, s, 2)), cluster[:,seqs[q][2]], [0.0, 0.0])*CompositeIndex(Index(seqs[q][3],FID{:f}(n, s, 1)), cluster[:,seqs[q][3]], [0.0, 0.0])
                ops = push!(ops, op, op')
            end
        end
    end
    ops_dict = Dict{Tuple{Int,Int}, Vector{Operator}}()
    hopping = Dict{Tuple{Int,Int}, Vector{Operator}}()
    getsites(op::Operator)=(op.id[1].index.site,op.id[2].index.site)
    for op in ops
        key = getsites(op)
        if !haskey(ops_dict, key)
            ops_dict[key] = []
        end
        push!(ops_dict[key], op)
    end
    for i = 1:length(cluster), j = 1:length(cluster)
        key = (i, j)
        hopping[key] = []
    end
    for key in keys(ops_dict)
        if haskey(hopping, key)
            hopping[key] = ops_dict[key]
        end
    end
    return hopping
end

"""
obtain the inverse cluster Green function with specific frequence ω : G⁻¹(ω)
"""
function invCGF(cluster::Cluster, bases::BinaryBases, hilbert::Hilbert, neighbors::Neighbors, ω::Number)
    

"""
obtain the CPT Green function with a periodic scheme : G    
"""




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
