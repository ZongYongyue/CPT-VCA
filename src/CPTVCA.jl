module CPTVCA
using Base.Threads
using Arpack: eigs
using LinearAlgebra: Eigen, dot, I
using SparseArrays: SparseMatrixCSC, spzeros, spdiagm
using QuantumLattices
using QuantumLattices: expand
using ExactDiagonalization
using ExactDiagonalization: eigen
using BlockArrays
using IterativeSolvers: minres, gmres
using StaticArrays: SVector

export Cluster, Clusters, CPTinfo, interclusterbonds, interhoppingm, CGF, CPTGF, CPTspec

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
The constant information needed in caculating the CPT green function
"""
struct CPTinfo
    cluster::Cluster
    clusters::Vector{Cluster}
    terms::Tuple{Vararg{Term}}
    hilbert::Hilbert
    neighbors::Neighbors
    spin::Rational 
    norbital::Int 
    table::Table 
    bases::BinaryBases
    ed_Hilbert::ED
    eigensystem::Eigen
    Ev::Matrix
    bases₁::BinaryBases
    bases₂::BinaryBases
    ed₁::ED
    ed₂::ED
    Hm₁::SparseMatrixCSC
    Hm₂::SparseMatrixCSC
    I₁::SparseMatrixCSC
    I₂::SparseMatrixCSC
    E₁::SparseMatrixCSC
    E₂::SparseMatrixCSC
    unitcellvectors::AbstractVector{<:AbstractVector{<:Number}}
    BZsteps::Tuple{Int, Vararg{Int}}
    η::Number 
end
function CPTinfo(lattice::AbstractLattice,bases::BinaryBases,terms::Tuple{Vararg{Term}},hilbert::Hilbert,neighbors::Neighbors,supervectors::AbstractVector{<:AbstractVector{<:Number}},unitcellvectors::AbstractVector{<:AbstractVector{<:Number}},BZsteps::Tuple{Int, Vararg{Int}},zeroplus::Number=0.05,steps::Tuple{Int, Vararg{Int}}=(1,1))
    cluster = Cluster(lattice)
    clusters = Clusters(cluster,steps,supervectors)
    spin = Rational((hilbert[1].nspin-1)/2)
    norbital = hilbert[1].norbital
    table = Table(hilbert, Metric(EDKind{:FED}(),hilbert))
    ed_Hilbert = ED(lattice, hilbert, terms, TargetSpace(bases))
    eigensystem = eigen(matrix(ed_Hilbert); nev=1)
    Ev = eigensystem.vectors
    bases₁=BinaryBases(hilbert[1].nspin*norbital*length(lattice),Int(bases.id[1][2]+1))
    bases₂=BinaryBases(hilbert[1].nspin*norbital*length(lattice),Int(bases.id[1][2]-1))
    ed₁ = ED(lattice, hilbert, terms, TargetSpace(bases₁))
    Hm₁ = expand(ed₁.Hₘ).contents[(TargetSpace(bases₁)[1],TargetSpace(bases₁)[1])].matrix
    ed₂ = ED(lattice, hilbert, terms, TargetSpace(bases₂))
    Hm₂ = expand(ed₂.Hₘ).contents[(TargetSpace(bases₂)[1],TargetSpace(bases₂)[1])].matrix
    I₁ = SparseMatrixCSC(spdiagm(0 => ones(length(bases₁))))
    I₂ = SparseMatrixCSC(spdiagm(0 => ones(length(bases₂))))
    E₁ = eigensystem.values[1]*I₁
    E₂ = eigensystem.values[1]*I₂ 
    η = zeroplus*im
    return CPTinfo(cluster,clusters,terms,hilbert,neighbors,spin,norbital,table,bases,ed_Hilbert,eigensystem,Ev,bases₁,bases₂,ed₁,ed₂,Hm₁,Hm₂,I₁,I₂,E₁,E₂,unitcellvectors,BZsteps,η)
end


"""
obtain the inter-cluster bonds between the given cluster and its surrounding clusters 
"""  
function interclusterbonds(info::CPTinfo)
    cluster = info.cluster
    neighbors = info.neighbors
    surroundings = filter(x->x!=cluster,info.clusters)
    for (index₁, val) in enumerate(info.clusters)
        if val == cluster
            cluster = val
        end
    end
    interclbs_dict =  Dict{Tuple,Bond}() 
    for (index₂, surrcluster) in enumerate(surroundings)
        seqs = interlinks(surrcluster, cluster, neighbors)
        if !isempty(seqs)
            for q in 1:length(seqs)
                key = ((cluster.id, seqs[q][2]),(surrcluster.id,seqs[q][3]))#((0,a),(m,b))
                Point₁ = Point(seqs[q][2],cluster[:,seqs[q][2]],cluster[:,1])
                Point₂ = Point(seqs[q][3],surrcluster[:,seqs[q][3]],surrcluster[:,1])
                interclbs_dict[key] = Bond(seqs[q][1], Point₁, Point₂)
            end
        end
    end
    return interclbs_dict
end 

"""
obtain the inter-cluster hopping matrix with specific momentum k
"""
function interhoppingm(info::CPTinfo, k::SVector)
    cluster = info.cluster
    terms = info.terms
    hilbert = info.hilbert
    interclbs_dict = interclusterbonds(info)
    spin = info.spin
    norbital = info.norbital
    hpv_dict = Dict{Tuple,Number}()
    getsites(op::Operator)=(Int(op.id[1].index.iid.spin+spin+1), Int(op.id[2].index.iid.spin+spin+1), op.id[1].index.iid.orbital, op.id[2].index.iid.orbital, op.id[1].index.site, op.id[2].index.site)
    getvalue(op::Operator)=op.value
    for (key₁, bond) in interclbs_dict
        ops = expand(terms[1], bond, hilbert, half=true)
        for (key₂, op) in ops.contents
            key₃ = getsites(op)
            hpv_dict[key₃] = getvalue(op)*exp(im*dot(k, op.id[2].icoordinate))
        end
    end
    A = zeros(Complex,Int(2*spin)+1,Int(2*spin)+1,norbital,norbital,length(cluster),length(cluster))
    B = zeros(Complex,(Int(2*spin)+1)*norbital*length(cluster),(Int(2*spin)+1)*norbital*length(cluster))
    C = BlockArray(B, fill(Int(2*spin)+1, norbital*length(cluster)), fill(Int(2*spin)+1, norbital*length(cluster)))
    for k in 1:norbital, l in 1:norbital, m in 1:length(cluster), n in 1:length(cluster)
        for i in 1:Int(2*spin)+1, j in 1:Int(2*spin)+1
            if (i,j,k,l,m,n) in keys(hpv_dict)
                A[i,j,k,l,m,n] = hpv_dict[(i,j,k,l,m,n)]
            else
                A[i,j,k,l,m,n] = Complex(0.0)
            end
            
        end
        C[Block((m-1)*norbital+k,(n-1)*norbital+l)]=A[:,:,k,l,m,n]
    end
    return Matrix(reshape(C,(Int(2*spin)+1)*norbital*length(cluster), (Int(2*spin)+1)*norbital*length(cluster)))
end

"""
obtain the cluster green function of sepcific frequence ω
"""
function CGF(info::CPTinfo, ω::Number)
    cluster = info.cluster
    terms = info.terms
    spin = info.spin
    norbital = info.norbital
    table = info.table
    Ev = info.Ev
    bases = info.bases
    bases₁ = info.bases₁
    bases₂ = info.bases₂
    eCGF_den = (ω + info.η + terms[2].value/2)*info.I₁ - info.Hm₁ + info.E₁
    hCGF_den = (ω + info.η + terms[2].value/2)*info.I₂ + info.Hm₂ - info.E₂
    cgf_dict = Dict{Tuple, Number}()
    for a in 1:length(cluster), b in  1:length(cluster)
        for  s1 in -spin:spin, s2 in -spin:spin, n1 in 1:norbital, n2 in 1:norbital
            key = (Int(s1+spin+1),Int(s2+spin+1),n1,n2,a,b)
            anniop = 1*CompositeIndex(Index(key[5], FID{:f}(n1, s1, 1)), cluster[:,key[5]], [0.0, 0.0])
            creatop = 1*CompositeIndex(Index(key[6], FID{:f}(n2, s2, 2)), cluster[:,key[6]], [0.0, 0.0])
            anniopm_left = matrix(anniop, (bases, bases₁), table)
            creatopm_right = matrix(creatop, (bases₁, bases), table)
            anniopm_right = matrix(anniop, (bases₂, bases), table)
            creatopm_left = matrix(creatop, (bases, bases₂), table)
            ecgf = dot((conj(transpose(Ev))*anniopm_left),gmres(eCGF_den,(creatopm_right*Ev))) 
            hcgf = dot((conj(transpose(Ev))*creatopm_left),gmres(hCGF_den,(anniopm_right*Ev)))
            cgf_dict[key] = ecgf + hcgf
        end
    end
    A = zeros(Complex,Int(2*spin)+1,Int(2*spin)+1,norbital,norbital,length(cluster),length(cluster))
    B = zeros(Complex,(Int(2*spin)+1)*norbital*length(cluster),(Int(2*spin)+1)*norbital*length(cluster))
    C = BlockArray(B, fill(Int(2*spin)+1, norbital*length(cluster)), fill(Int(2*spin)+1, norbital*length(cluster)))
    for k in 1:norbital, l in 1:norbital, m in 1:length(cluster), n in 1:length(cluster)
        for i in 1:Int(2*spin)+1, j in 1:Int(2*spin)+1
            if (i,j,k,l,m,n) in keys(cgf_dict)
                A[i,j,k,l,m,n] = cgf_dict[(i,j,k,l,m,n)]
            else
                A[i,j,k,l,m,n] = Complex(0.0)
            end
            
        end
        C[Block((m-1)*norbital+k,(n-1)*norbital+l)]=A[:,:,k,l,m,n]
    end
    return Matrix(reshape(C,(Int(2*spin)+1)*norbital*length(cluster), (Int(2*spin)+1)*norbital*length(cluster))) 
end

"""
obtaion the CPT green function with (k, ω)
"""
function CPTGF(info::CPTinfo, k::SVector, ω::Number)
    cluster = info.cluster
    spin = info.spin
    norbital = info.norbital
    V = interhoppingm(info, k)
    G = CGF(info, ω)
    In = Matrix{Complex}(I, size(G,1), size(G,2))
    GFm = G*inv(In-V*G)
    GFB = BlockArray(GFm, fill((Int(2*spin)+1)*norbital, length(cluster)),fill((Int(2*spin)+1)*norbital, length(cluster)))
    GF = 0
    for i in 1:length(cluster), j in 1:length(cluster)
        GF += sum(GFB[Block(i,j)])*exp(-im*dot(k, (cluster[:,i] -cluster[:,j])))
    end
    return (1/length(cluster))*GF
end

"""
obtain the spectral function of CPT green function 
"""
@inline function CPTspec(info, k_path, ω_range)
    path = selectpath(BrillouinZone(reciprocals(info.unitcellvectors), info.BZsteps), k_path)
    A = zeros(Float64, length(ω_range), length(path[1]))
    function calculate_element(m,i)
        k = path[1][i]
        ω = ω_range[m]
        return (-1 / π) * imag(CPTGF(info, k, ω))
    end
    @threads for i in eachindex(path[1])
        for m in eachindex(ω_range)
            A[m, i] = calculate_element(m, i)
        end
    end
    return A
end


end # module CPTVCA
