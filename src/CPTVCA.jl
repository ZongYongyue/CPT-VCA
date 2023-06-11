module CPTVCA
using Base.Threads
using Arpack: eigs
using LinearAlgebra: Eigen, dot, I, SymTridiagonal, Hermitian
using SparseArrays: SparseMatrixCSC, spzeros, spdiagm
using QuantumLattices
using QuantumLattices: expand
using ExactDiagonalization
using ExactDiagonalization: eigen
using BlockArrays
using IterativeSolvers: minres, gmres
using StaticArrays: SVector
using KrylovKit
export Cluster, Clusters, CPTinfo, gen_krylov, interclusterbonds, interhoppingm, CGF, CPTGF, CPTspec

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
    hoppingterms::Tuple{Vararg{Term}}
    Hubbardterm::Term
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
function CPTinfo(lattice::AbstractLattice,bases::BinaryBases,terms::Tuple{Vararg{Term}},hoppingterms::Tuple{Vararg{Term}}, Hubbardterm::Term,hilbert::Hilbert,neighbors::Neighbors,supervectors::AbstractVector{<:AbstractVector{<:Number}},unitcellvectors::AbstractVector{<:AbstractVector{<:Number}},BZsteps::Tuple{Int, Vararg{Int}},zeroplus::Number=0.05,steps::Tuple{Int, Vararg{Int}}=(1,1))
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
    return CPTinfo(cluster,clusters,terms,hoppingterms,Hubbardterm,hilbert,neighbors,spin,norbital,table,bases,ed_Hilbert,eigensystem,Ev,bases₁,bases₂,ed₁,ed₂,Hm₁,Hm₂,I₁,I₂,E₁,E₂,unitcellvectors,BZsteps,η)
end

"""
calculate the Lanczos bases and tridiagonal matrix in krylov subspace
"""
function gen_krylov(Matrix, invec, m)
    Matrix = Hermitian(Matrix)
    orth = KrylovKit.ModifiedGramSchmidt()
    iterator = LanczosIterator(Matrix, invec, orth)
    factorization = KrylovKit.initialize(iterator)
    for _ in 1:m-1
        KrylovKit.expand!(iterator, factorization)
    end
    basis_vectors = basis(factorization)
    T = rayleighquotient(factorization)
    return basis_vectors, T
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
calculate the inter-cluster hopping matrix with specific momentum k
"""
function interhoppingm(info::CPTinfo, k::SVector)
    hoppingterms = info.hoppingterms
    cluster = info.cluster
    hilbert = info.hilbert
    interclbs_dict = interclusterbonds(info)
    spin = info.spin
    norbital = info.norbital
    hpv_dict = Dict{Tuple,Number}()
    getsites(op::Operator)=(Int(op.id[1].index.iid.spin+spin+1), Int(op.id[2].index.iid.spin+spin+1), op.id[1].index.iid.orbital, op.id[2].index.iid.orbital, op.id[1].index.site, op.id[2].index.site)
    getvalue(op::Operator)=op.value
    for (key₁, bond) in interclbs_dict
        for t in 1:length(hoppingterms)
            ops = expand(hoppingterms[t], bond, hilbert, half=true)
            for (key₂, op) in ops.contents
                key₃ = getsites(op)
                hpv_dict[key₃] = getvalue(op)*exp(im*dot(k, op.id[2].icoordinate))
            end
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
calculate the cluster green function of sepcific frequence ω
"""
function CGF(info::CPTinfo, ω::Number)
    Hubbardterm = info.Hubbardterm
    cluster = info.cluster
    spin = info.spin
    norbital = info.norbital
    table = info.table
    Ev = info.Ev
    bases = info.bases
    bases₁ = info.bases₁
    bases₂ = info.bases₂
    eCGF_den = (ω + info.η + Hubbardterm.value/2)*info.I₁ - info.Hm₁ + info.E₁
    hCGF_den = (ω + info.η + Hubbardterm.value/2)*info.I₂ + info.Hm₂ - info.E₂
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
obtain the cluster green function with Lanczos method 
"""
function CGF(info::CPTinfo, ω_range::StepRangeLen)
    Hubbardterm = info.Hubbardterm
    cluster = info.cluster
    spin = info.spin
    norbital = info.norbital
    table = info.table
    Ev = info.Ev
    bases = info.bases
    bases₁ = info.bases₁
    bases₂ = info.bases₂
    In = Matrix{Complex}(I, 300, 300)
    ecgf_dict = Dict{Tuple, Vector}()
    hcgf_dict = Dict{Tuple, Vector}()
    Threads.@threads for index₁ in 1:length(cluster)
        for s₁ in -spin:spin, n₁ in 1:norbital
            Cr⁺ = 1*CompositeIndex(Index(index₁, FID{:f}(n₁, s₁, 2)), cluster[:, index₁], [0.0, 0.0])
            Cr  = 1*CompositeIndex(Index(index₁, FID{:f}(n₁, s₁, 1)), cluster[:, index₁], [0.0, 0.0])
            Cr⁺Ω = (matrix(Cr⁺, (bases₁, bases), table)*Ev)[:,1]
            CrΩ = (matrix(Cr, (bases₂, bases), table)*Ev)[:,1]
            krylovbasis_e, T_e = gen_krylov(info.Hm₁, Cr⁺Ω, 300)
            krylovbasis_h, T_h = gen_krylov(info.Hm₂, CrΩ, 300)
            for index₂ in 1:length(cluster)
                for s₂ in -spin:spin, n₂ in 1:norbital
                    key₁ = (Int(s₂+spin+1),Int(s₁+spin+1),n₂,n₁,index₂,index₁)
                    key₂ = (Int(s₁+spin+1),Int(s₂+spin+1),n₁,n₂,index₁,index₂)
                    Cl⁺ = 1*CompositeIndex(Index(index₂, FID{:f}(n₂, s₂, 2)), cluster[:, index₂], [0.0, 0.0])
                    Cl  = 1*CompositeIndex(Index(index₂, FID{:f}(n₂, s₂, 1)), cluster[:, index₂], [0.0, 0.0])
                    Cl⁺Ω = (matrix(Cl⁺, (bases₁, bases), table)*Ev)[:,1]
                    ClΩ = (matrix(Cl, (bases₂, bases), table)*Ev)[:,1]
                    X_e = KrylovKit.project!(zeros(Complex, 300), krylovbasis_e, Cl⁺Ω)
                    X_h = KrylovKit.project!(zeros(Complex, 300), krylovbasis_h, ClΩ)
                    egf_vec = Vector()
                    hgf_vec = Vector()
                    for ω in ω_range
                        egf = dot(X_e, inv((ω + info.η + Hubbardterm.value/2 + info.eigensystem.values[1])*In - Matrix(T_e))[:,1])*√(Cr⁺Ω'*Cr⁺Ω)
                        hgf = dot(X_h, inv((ω + info.η + Hubbardterm.value/2 - info.eigensystem.values[1])*In + Matrix(T_h))[:,1])*√(CrΩ'*CrΩ)
                        push!(egf_vec, egf)
                        push!(hgf_vec, hgf)
                    end
                    ecgf_dict[key₁] = egf_vec
                    hcgf_dict[key₂] = hgf_vec
                end
            end
        end
    end
    cgf_vec = Vector{Matrix}()
    for ω in eachindex(ω_range)
        A = zeros(Complex,Int(2*spin)+1,Int(2*spin)+1,norbital,norbital,length(cluster),length(cluster))
        B = zeros(Complex,(Int(2*spin)+1)*norbital*length(cluster),(Int(2*spin)+1)*norbital*length(cluster))
        C = BlockArray(B, fill(Int(2*spin)+1, norbital*length(cluster)), fill(Int(2*spin)+1, norbital*length(cluster)))
        for k in 1:norbital, l in 1:norbital, m in 1:length(cluster), n in 1:length(cluster)
            for i in 1:Int(2*spin)+1, j in 1:Int(2*spin)+1
                if (i,j,k,l,m,n) in keys(ecgf_dict)
                    A[i,j,k,l,m,n] = ecgf_dict[(i,j,k,l,m,n)][ω] + hcgf_dict[(i,j,k,l,m,n)][ω]
                else
                    A[i,j,k,l,m,n] = Complex(0.0)
                end
                
            end
            C[Block((m-1)*norbital+k,(n-1)*norbital+l)]=A[:,:,k,l,m,n]
        end
        push!(cgf_vec, Matrix(reshape(C,(Int(2*spin)+1)*norbital*length(cluster), (Int(2*spin)+1)*norbital*length(cluster)))) 
    end
    return cgf_vec
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
obtaion the CPT green function with Lanczos Method 
"""
function CPTGF(info::CPTinfo, k_path::Tuple, ω_range::StepRangeLen)
    path = selectpath(BrillouinZone(reciprocals(info.unitcellvectors), info.BZsteps), k_path)
    G = zeros(Complex, length(ω_range), length(path[1]))
    cluster = info.cluster
    spin = info.spin
    norbital = info.norbital
    V_vec = Vector{Matrix}()
    for k in path[1]
        push!(V_vec, interhoppingm(info, k))
    end
    G_vec = CGF(info, ω_range)
    In = Matrix{Complex}(I, size(G_vec[1],1), size(G_vec[1],2))
    for i in eachindex(ω_range) , j in eachindex(path[1])
        k = path[1][j]
        GFm = G_vec[i]*inv(In-V_vec[j]*G_vec[i])
        GFB = BlockArray(GFm, fill((Int(2*spin)+1)*norbital, length(cluster)),fill((Int(2*spin)+1)*norbital, length(cluster)))
        GF = 0
        for m in 1:length(cluster), n in 1:length(cluster)
            GF += sum(GFB[Block(m,n)])*exp(-im*dot(k, (cluster[:,m] -cluster[:,n])))
        end
        G[i, j] = (1/length(cluster))*GF
    end
    return G
end


"""
obtain the spectral function of CPT green function 
"""
@inline function CPTspec(info, k_path, ω_range; lanczos::Bool=false)
    if lanczos
        G = CPTGF(info, k_path, ω_range)
        A = (-1 / π) * imag(G)
    else
        path = selectpath(BrillouinZone(reciprocals(info.unitcellvectors), info.BZsteps), k_path)
        A = zeros(Float64, length(ω_range), length(path[1]))
        function calculate_element(m, i)
            k = path[1][i]
            ω = ω_range[m]
            return (-1 / π) * imag(CPTGF(info, k, ω))
        end
        Threads.@threads for i in eachindex(path[1])
            for m in eachindex(ω_range)
                A[m, i] = calculate_element(m, i)
            end
        end
    end
    return A
end




end # module CPTVCA