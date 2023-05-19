module CPTVCA

using Arpack: eigs
using LinearAlgebra: Eigen, dot
using SparseArrays: SparseMatrixCSC, spzeros, spdiagm
using QuantumLattices
using ExactDiagonalization
using ExactDiagonalization: eigen
using BlockArrays
using IterativeSolvers: minres

export Cluster, Clusters, interhopping, ClusterH, CGF, CPTGF, specGF

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
obtaion a supercluster consisting of a given cluster and its surrounding clusters 
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
construct a type indexed by cluster sites (a, b) that contain the information of operators, the matrix of the operators and the value of given state
"""
struct HPMatrixRepresentation
    ops::OperatorSum
    opsmatrix::SparseMatrixCSC{ComplexF64, Int64}
    value::Number
end 

"""
obtain the inter-cluster hopping  between the given cluster and its surrounding clusters with specific momentum k
"""
function interhopping(clusters::Vector{Cluster}, cluster::Cluster, eigensystem::Eigen, bases::BinaryBases, hilbert::Hilbert, neighbors::Neighbors, k::Vector)
    surroundings = filter(x->x!=cluster,clusters)
    spin = Rational((hilbert[1].nspin-1)/2)
    norbital = hilbert[1].norbital
    ops = Vector{Operator}()
    for (index,value) in enumerate(surroundings), s in -spin:spin, n in 1:norbital
        seqs = interlinks(value, cluster, neighbors)
        if !isempty(seqs)
            for q in 1:length(seqs)
                op = exp(-im*dot(k, value[:,1]))*CompositeIndex(Index(seqs[q][2],FID{:f}(n, s, 2)), cluster[:,seqs[q][2]], [0.0, 0.0])*CompositeIndex(Index(seqs[q][3],FID{:f}(n, s, 1)), cluster[:,seqs[q][3]], [0.0, 0.0])
                ops = push!(ops, op)
            end
        end
    end
    ops_dict = Dict{Tuple{Int,Int}, Vector{Operator}}()
    hpm_dict = Dict{Tuple{Int,Int}, HPMatrixRepresentation}() 
    hpv_dict = Dict{Tuple{Int,Int}, Number}()
    getsites(op::Operator)=(op.id[1].index.site,op.id[2].index.site)
    for op in ops
        key = getsites(op)
        if !haskey(ops_dict, key)
            ops_dict[key] = []
        end
        push!(ops_dict[key], op)
    end
    for key in keys(ops_dict)
        ops_hasconju = sum(ops_dict[key])+sum(ops_dict[key])'
        hoppingmatrix = matrix(ops_hasconju, (bases,bases), Table(hilbert, Metric(EDKind{:FED}(),hilbert)))
        hpvalue = dot(eigensystem.vectors,hoppingmatrix*eigensystem.vectors)
        #hpm_dict[key] = HPMatrixRepresentation(ops_hasconju, hoppingmatrix, hpvalue) # can choose to return this 
        hpv_dict[key] = hpvalue
    end
    return hpv_dict
end

"""
obtaion the cluster Green function with specific frequence
"""
function CGF(cluster::Cluster, bases_Hilbert::BinaryBases, bases_Fock::BinaryBases, eigensystem_Hilbert::Eigen, Hmatrix_Fock::SparseMatrixCSC, hilbert::Hilbert, ω::Number)
    table = Table(hilbert, Metric(EDKind{:FED}(),hilbert))
    Iₘ = SparseMatrixCSC(spdiagm(0 => ones(size(Hmatrix_Fock)[1])))
    Evalue = eigensystem_Hilbert.values[1]*Iₘ
    ωₘ = ω*Iₘ 
    Evector = eigensystem_Hilbert.vectors
    eCGF_den = ωₘ - Hmatrix_Fock + Evalue
    hCGF_den = ωₘ + Hmatrix_Fock - Evalue
    spin = Rational((hilbert[1].nspin-1)/2)
    norbital = hilbert[1].norbital
    cgf_dict = Dict{Tuple{Int,Int}, Number}()
    for a in 1:length(cluster), b in  1:length(cluster)
        key = (a, b)
        cgf_dict[key] = 0
        for  s in -spin:spin, n in 1:norbital
            anniop = 1*CompositeIndex(Index(key[1], FID{:f}(n, s, 1)), cluster[:,key[1]], [0.0, 0.0])
            creatop = 1*CompositeIndex(Index(key[2], FID{:f}(n, s, 2)), cluster[:,key[2]], [0.0, 0.0])
            anniopm_left = matrix(anniop, (bases_Hilbert, bases_Fock), table)
            creatopm_right = matrix(creatop, (bases_Fock, bases_Hilbert), table)
            anniopm_right = matrix(anniop, (bases_Fock, bases_Hilbert), table)
            creatopm_left = matrix(creatop, (bases_Hilbert, bases_Fock), table)
            cgf_dict[key] += dot(conj(transpose(Evector))*anniopm_left,minres(eCGF_den,creatopm_right*Evector)) + dot(conj(transpose(Evector))*creatopm_left,minres(hCGF_den,anniopm_right*Evector))
        end
    end
    return cgf_dict
end

"""
obtain the matrix representation of Hamiltonian for a single cluster
"""
function ClusterH(lattice::AbstractLattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, targetspace::TargetSpace; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)
    ed = ED(lattice, hilbert, terms, targetspace; neighbors, boundary)
    ClusterHₘ = expand(ed.Hₘ).contents[(targetspace[1],targetspace[1])].matrix
    return ClusterHₘ
end 

"""
obtain the CPT Green function with a periodic scheme   
"""
function CPTGF(lattice::AbstractLattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, bases::BinaryBases, neighbors::Neighbors, k::Vector, ω::Number, vectors::AbstractVector{<:AbstractVector{<:Number}}, steps::Tuple{Int, Vararg{Int}}=(1,1))
    cluster = Cluster(lattice)
    clusters = Clusters(cluster, steps, vectors)
    ed_Hilbert = ED(lattice, hilbert, terms, TargetSpace(bases))
    eigensystem = eigen(matrix(ed_Hilbert); nev=1)
    nstates = 0
    for (index, values) in enumerate(hilbert)
        nstates += hilbert[index].norbital*hilbert[index].nspin
    end
    bases_f = BinaryBases(nstates)
    Hm_f = ClusterH(lattice, hilbert, terms, TargetSpace(bases_f))
    hoppingvalues =  interhopping(clusters, cluster, eigensystem, bases, hilbert, neighbors, k)
    cgfvalues = CGF(cluster, bases, bases_f, eigensystem, Hm_f, hilbert, ω)
    cptgf_dict = Dict{Tuple{Int,Int}, Number}()
    for (key, value) in cgfvalues
        if haskey(hoppingvalues, key)
            cptgf_dict[key] = 1/value - hoppingvalues[key]
        else
            cptgf_dict[key] = 1/value
        end
    end
    cptgf = 0
    for (key, value) in cptgf_dict
        cptgf += exp(-im*dot(k, (cluster[:,key[1]] -cluster[:,key[2]])))*(1/value)
    end
    return (1/length(lattice))*cptgf
end

end # module CPTVCA
