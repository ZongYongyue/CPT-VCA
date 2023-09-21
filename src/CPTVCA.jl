module CPTVCA

using Arpack
using Distributed
using ExactDiagonalization
using ExactDiagonalization: matrix
using KrylovKit
using LinearAlgebra
using QuadGK
using QuantumLattices
using Serialization
using TimerOutputs

export Conservation, VCA, VCAs, selectBases, clusterGreenFunction, clusterGreenFunctionLoop, origiQuadraticTerms!, referQuadraticTerms, EDSolver, GPintegrand, GPcore, origiQuadraticTerms
export perGreenFunction, GreenFunctionPath, singleParticleGreenFunction, clusterspectrum, spectrum, Optimal, optimals, OrderParameters, GrandPotential, spawn, saveData, loadData
export antiferro, belongs, moirephase, spinrotation, Perioder
export initialstate, showgpint, showvm

const vcatimer = TimerOutput()

"""
    abstract type InitialState

The Operators needed in calculating the initial vector in Lanczos iteration method of krylov subspace
"""
abstract type InitialState end

"""
    FInitialState{O<:Operator} <: InitialState
    SInitialState{O<:Operator} <: InitialState

The Operators in a fermion/spin system needed in calculating the initial vector in Lanczos iteration method of krylov subspace 
"""
struct FInitialState{O<:Operator} <: InitialState
    sign::Int
    ops::OperatorSum{O}
    function FInitialState(sign::Int, key::Tuple)
        ops = Operators(1*CompositeIndex(Index(key[2], FID{:f}(key[3], key[1], sign)), [0.0, 0.0], [0.0, 0.0]))
        new{eltype(ops)}(sign, ops)
    end
end
struct SInitialState <: InitialState
    sign::Int
    key::Tuple
    function SInitialState(sign::Int, key::Tuple)
        #to be update
    end
end

"""
    initialstate(::EDKind{:FED}, sign::Int, key::Tuple) -> FInitialState
    initialstate(::EDKind{:SED}, sign::Int, key::Tuple) -> SInitialState

Get the Operators needed in calculating the initial vector
"""
function initialstate(::EDKind{:FED}, sign::Int, key::Tuple) FInitialState(sign, key) end
function initialstate(::EDKind{:SED}, sign::Int, key::Tuple) SInitialState(sign, key) end

"""
    Conservation{S<:Int, P<:Rational, O<:Union{Int,Nothing}, R<:Union{Real,Nothing}}

The information of conserved quantities
"""
struct Conservation{S<:Int, P<:Rational, O<:Union{Int,Nothing}, R<:Union{Real,Nothing}}
    nsite::S
    nspin::P
    nparticle::O
    spinz::R
    function Conservation(ns::Int, nspin::Rational; np::Union{Int,Nothing}=ns, sz::Union{Real,Nothing}=nothing)
        return new{typeof(ns),typeof(nspin),typeof(np),typeof(sz)}(ns, nspin, np, sz)
    end
end

"""
    selectBases(cons::Conservation)

Select a method to express bases considering the conserved quantities
"""
function selectBases(cons::Conservation)
    if isnothing(cons.nparticle)
        return BinaryBases(2cons.nsite)
    else
        if isnothing(cons.spinz)
            return BinaryBases(2cons.nsite, cons.nparticle)
        else
            return BinaryBases(1:cons.nsite, Int((cons.nparticle-2cons.spinz)/2))⊗BinaryBases(cons.nsite+1:2cons.nsite, Int((cons.nparticle+2cons.spinz)/2))
        end
    end
end

"""
    genkrylov(matrix::AbstractMatrix, sstate::AbstractVector, m::Int)

Generate a krylov subspace with Lanczos iteration method
"""
function genkrylov(matrix::AbstractMatrix, sstate::AbstractVector, m::Int)
    orth = KrylovKit.ModifiedGramSchmidt()
    iterator = LanczosIterator(matrix, sstate, orth)
    factorization = KrylovKit.initialize(iterator)
    for _ in 1:m-1
        KrylovKit.expand!(iterator, factorization)
    end
    basis_vectors = basis(factorization)
    T = rayleighquotient(factorization)
    a, c = [0.0; T.ev[1:end-1]], [T.ev[1:end-1]; 0.0]
    return basis_vectors, a, T.dv, c
end 

"""
    Kryvals{T<:AbstractVector, S<:AbstractVector, R<:AbstractVector, N<:AbstractVector, P<:AbstractMatrix}

The information obtained with the krylov subspace method that needed in calculating the cluster Green function
"""
struct Kryvals{T<:AbstractVector, S<:AbstractVector, R<:AbstractVector, N<:AbstractVector, P<:AbstractMatrix}
    avs::T
    bvs::S
    cvs::R
    norms::N
    projects::P
    function Kryvals(H::AbstractVector, initialsets::AbstractVector, m::Int=200)
        n = length(initialsets)
        avs, bvs, cvs, norms, projects = Vector{Vector{Float64}}(undef,n), Vector{Vector{Float64}}(undef,n), Vector{Vector{Float64}}(undef,n), Vector{Float64}(undef,n), Matrix{Vector{ComplexF64}}(undef,n,n)
        if length(H) == 2
            for i in 1:Int(n/2)
                krybasis, avs[i], bvs[i], cvs[i] = genkrylov(H[1], initialsets[i], m)  
                norms[i] = norm(initialsets[i])
                for j in 1:Int(n/2)
                    projects[i, j] = KrylovKit.project!(zeros(ComplexF64, m), krybasis, initialsets[j])
                end
                for j in Int(n/2)+1:n
                    projects[i, j] = zeros(ComplexF64, m)
                end
            end
            for i in Int(n/2)+1:n
                krybasis, avs[i], bvs[i], cvs[i] = genkrylov(H[2], initialsets[i], m)  
                norms[i] = norm(initialsets[i])
                for j in Int(n/2)+1:n
                    projects[i, j] = KrylovKit.project!(zeros(ComplexF64, m), krybasis, initialsets[j])
                end
                for j in 1:Int(n/2)
                    projects[i, j] = zeros(ComplexF64, m)
                end
            end
        elseif length(H) == 1
            for i in 1:n
                krybasis, avs[i], bvs[i], cvs[i] = genkrylov(H[1], initialsets[i], m)  
                norms[i] = norm(initialsets[i])
                for j in 1:n
                    projects[i, j] = KrylovKit.project!(zeros(ComplexF64, m), krybasis, initialsets[j])
                end
            end
        end
        new{typeof(avs), typeof(bvs), typeof(cvs), typeof(norms), typeof(projects)}(avs, bvs, cvs, norms, projects)
    end
end

"""
    initialsets(k::EDKind, gs::AbstractVector, cons::Conservation, table::Table)

Give two sets of initial states used in generating Krylov subspaces for advanced and retarded Green function
"""
function initialsets(k::EDKind, gs::AbstractVector, cons::Conservation, table::Table)
    sets₁, sets₂ = Vector{Vector}(undef, length(table)), Vector{Vector}(undef, length(table))
    orderkeys = sort(collect(keys(table)), by = x -> table[x])
    for i in eachindex(orderkeys)
        ops₁, ops₂ = initialstate(k, 1, orderkeys[i]).ops, initialstate(k, 2, orderkeys[i]).ops
        if !isnothing(cons.nparticle)&&!isnothing(cons.spinz)
            spin₁, spin₂ = collect(ops₁)[1].id[1].index.iid.spin, collect(ops₂)[1].id[1].index.iid.spin
            initialstate₁ = (matrix(ops₁, (selectBases(Conservation(cons.nsite, cons.nspin; np=cons.nparticle-1, sz=cons.spinz-spin₁)), selectBases(cons)), table)*gs)[:,1]
            initialstate₂ = (matrix(ops₂, (selectBases(Conservation(cons.nsite, cons.nspin; np=cons.nparticle+1, sz=cons.spinz+spin₂)), selectBases(cons)), table)*gs)[:,1]
        elseif !isnothing(cons.nparticle)&&isnothing(cons.spinz)
            initialstate₁ = (matrix(ops₁, (selectBases(Conservation(cons.nsite, cons.nspin; np=cons.nparticle-1)), selectBases(cons)), table)*gs)[:,1]
            initialstate₂ = (matrix(ops₂, (selectBases(Conservation(cons.nsite, cons.nspin; np=cons.nparticle+1)), selectBases(cons)), table)*gs)[:,1]
        else
            initialstate₁ = (matrix(ops₁, (selectBases(Conservation(cons.nsite, cons.nspin; np=nothing)), selectBases(cons)), table)*gs)[:,1]
            initialstate₂ = (matrix(ops₂, (selectBases(Conservation(cons.nsite, cons.nspin; np=nothing)), selectBases(cons)), table)*gs)[:,1]
        end
        sets₁[i] = initialstate₁ 
        sets₂[i] = initialstate₂
    end
    return sets₁, sets₂
end

"""
    Sysvals{K<:EDKind, R<:Number}

The all information needed to calculate the Green function of a finite size system 
"""
struct Sysvals{K<:EDKind, R<:Number}
    gsenergy::R
    kryvals₁::Kryvals
    kryvals₂::Kryvals
end
function Sysvals(k::EDKind, gse::Real, gs::AbstractVector, rops::OperatorSum, cons::Conservation, table::Table; m::Int=200)
    sets₁, sets₂ = initialsets(k, gs, cons, table)
    if !isnothing(cons.nparticle)&&!isnothing(cons.spinz)
        H₁1 = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=cons.nparticle-1,sz=cons.spinz+cons.nspin)),selectBases(Conservation(cons.nsite,cons.nspin;np=cons.nparticle-1,sz=cons.spinz+cons.nspin))), table)
        H₁2 = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=cons.nparticle-1,sz=cons.spinz-cons.nspin)),selectBases(Conservation(cons.nsite,cons.nspin;np=cons.nparticle-1,sz=cons.spinz-cons.nspin))), table)
        H₂1 = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=cons.nparticle+1,sz=cons.spinz-cons.nspin)),selectBases(Conservation(cons.nsite,cons.nspin;np=cons.nparticle+1,sz=cons.spinz-cons.nspin))), table)
        H₂2 = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=cons.nparticle+1,sz=cons.spinz+cons.nspin)),selectBases(Conservation(cons.nsite,cons.nspin;np=cons.nparticle+1,sz=cons.spinz+cons.nspin))), table)
        kryvals₁, kryvals₂ = Kryvals([H₁1,H₁2], sets₁, m), Kryvals([H₂1,H₂2], sets₂, m)
    elseif !isnothing(cons.nparticle)&&isnothing(cons.spinz)
        H₁ = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=cons.nparticle-1)),selectBases(Conservation(cons.nsite,cons.nspin;np=cons.nparticle-1))), table)
        H₂ = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=cons.nparticle+1)),selectBases(Conservation(cons.nsite,cons.nspin;np=cons.nparticle+1))), table)
        kryvals₁, kryvals₂ = Kryvals([H₁], sets₁, m), Kryvals([H₂], sets₂, m)
    else
        H₁ = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=nothing)),selectBases(Conservation(cons.nsite,cons.nspin;np=nothing))), table)
        H₂ = matrix(rops, (selectBases(Conservation(cons.nsite, cons.nspin;np=nothing)),selectBases(Conservation(cons.nsite,cons.nspin;np=nothing))), table)
        kryvals₁, kryvals₂ = Kryvals([H₁], sets₁, m), Kryvals([H₂], sets₂, m)
    end
    return Sysvals{typeof(k), typeof(gse)}(gse, kryvals₁, kryvals₂)
end

"""
    The solver to calculate the cluster Green function
"""
abstract type GFSolver end

"""
    struct EDSolver <: GFSolver 

The ED solver to calculate the cluster Green function
"""
struct EDSolver <: GFSolver 
    sysvals::Sysvals
end

function EDSolver(k::EDKind, refergenerator::OperatorGenerator, cons::Conservation, table::Table; m::Int=200)
    rops = expand(refergenerator)
    Hₘ = matrix(rops, (selectBases(cons), selectBases(cons)), table)
    vals, vecs, _  = eigsolve(Hₘ, 1, :SR, Float64)
    gse, gs = real(vals[1]), vecs[1]
    sysvals = Sysvals(k, gse, gs, rops, cons, table; m=m)
    EDSolver(sysvals)
end
"""
    thomas(a::AbstractVector, b::AbstractVector, c::AbstractVector, d::AbstractVector, n::Int)

Thomas method, used to solve linear equations respect to tridiagonal matrix
"""
function thomas(a::AbstractVector, b::AbstractVector, c::AbstractVector, d::AbstractVector, n::Int)
    x = Complex.(d)
    c_prime = Complex.(c)
    c_prime[1] /= b[1]
    x[1] /= b[1]
    for i = 2:n
        scale = 1.0 / (b[i] - c_prime[i-1]*a[i])
        c_prime[i] *= scale
        x[i] = (x[i] - a[i] * x[i-1]) * scale
    end
    for i = n-1:-1:1
        x[i] -= (c_prime[i] * x[i+1])
    end
    return x
end

"""
    solvercore(a::Complex, r::Complex, avs₁::AbstractVector, bvs₁::AbstractVector, cvs₁::AbstractVector, norm₁::AbstractVector, proj₁::AbstractMatrix, avs₂::AbstractVector, bvs₂::AbstractVector, cvs₂::AbstractVector, norm₂::AbstractVector, proj₂::AbstractMatrix, d::AbstractVector, m::Int, gfm::AbstractMatrix)

The calculating core of the method clusterGreenFunction
"""
function solvercore(a::Complex, r::Complex, avs₁::AbstractVector, bvs₁::AbstractVector, cvs₁::AbstractVector, norm₁::AbstractVector, proj₁::AbstractMatrix, avs₂::AbstractVector, bvs₂::AbstractVector, cvs₂::AbstractVector, norm₂::AbstractVector, proj₂::AbstractMatrix, d::AbstractVector, m::Int, gfm::AbstractMatrix)
    sgfm = copy(gfm)
    bv₁, bv₂ = [a .+ bvs₁[i] for i in 1:length(norm₁)], [r .- bvs₂[i] for i in 1:length(norm₂)]
    tmpv₁, tmpv₂ = [thomas(avs₁[i],bv₁[i],cvs₁[i],d,m) for i in 1:length(norm₁)], [thomas(-avs₂[i],bv₂[i],-cvs₂[i],d,m) for i in 1:length(norm₂)]
    for i in eachindex(norm₁)
        for j in eachindex(norm₂)
            sgfm[i, j] = dot(proj₁[i, j], tmpv₁[i]) * norm₁[i] + dot(proj₂[j, i], tmpv₂[j]) * norm₂[j]
        end
    end
    return sgfm
end
"""
    clusterGreenFunction(solver::GFSolver, ω::Complex)

Calculate the cluster Green function with krylov subspace method respect to single ω
"""
function clusterGreenFunction(solver::GFSolver, ω::Complex)
    kry₁, kry₂ = solver.sysvals.kryvals₁, solver.sysvals.kryvals₂
    avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂ = kry₁.avs, kry₁.bvs, kry₁.cvs, kry₁.norms, kry₁.projects, kry₂.avs, kry₂.bvs, kry₂.cvs, kry₂.norms, kry₂.projects
    gfm = zeros(ComplexF64, length(norm₁), length(norm₁))
    m = length(avs₁[1])
    d = [1.0;zeros(m-1)]
    a, r = ω - solver.sysvals.gsenergy, ω + solver.sysvals.gsenergy
    return solvercore(a, r, avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂, d, m, gfm)
end

"""
    clusterGreenFunctionLoop(solver::GFSolver, ω_range::AbstractRange)

Calculate the cluster Green function with krylov subspace method respect to a rnage of ω
"""
function clusterGreenFunctionLoop(solver::GFSolver, ω_range::AbstractRange)
    kry₁, kry₂ = solver.sysvals.kryvals₁, solver.sysvals.kryvals₂
    avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂ = kry₁.avs, kry₁.bvs, kry₁.cvs, kry₁.norms, kry₁.projects, kry₂.avs, kry₂.bvs, kry₂.cvs, kry₂.norms, kry₂.projects
    cgfv = [Matrix{ComplexF64}(undef, length(norm₁), length(norm₁)) for _ in 1:length(ω_range)]
    gfm = zeros(ComplexF64, length(norm₁), length(norm₁))
    m = length(avs₁[1])
    d = [1.0;zeros(m-1)]
    for i in eachindex(ω_range)
        a, r = ω_range[i] - solver.sysvals.gsenergy, ω_range[i] + solver.sysvals.gsenergy
        cgfv[i] = solvercore(a, r, avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂, d, m, gfm)
    end
    return cgfv
end

"""
以上可以算作ED部分的扩展
===============================================================
下面是和CPTVCA相关内容
"""

"""
    Perioder{P<:AbstractVector{<:Integer}, T<:AbstractArray{P}, S<:AbstractArray{P}, C<:AbstractVector{<:Tuple}}
    Perioder(unitcell::AbstractLattice, cluster::AbstractLattice, table::Table)

User should ensure that the cluster you choosed is compatible with the lattice generated by the unitcell you input and the unitcell you input should be enclosed in the cluster you choosed sharing an original point with the cluster.
"""
struct Perioder{P<:AbstractVector{<:Integer}, T<:AbstractArray{P}, S<:AbstractArray{P}, C<:AbstractVector{<:Tuple}}
    map₁::T
    map₂::S
    channels::C
end 
function Perioder(unitcell::AbstractLattice, cluster::AbstractLattice, table::Table)
    @assert !isempty(unitcell.vectors) "the vectors in unitcell cannot be empty !"
    seq = sort(collect(keys(table)), by = x -> table[x])
    nspin, norbi = sort(collect(Set(key[1] for key in seq))), sort(collect(Set(key[3] for key in seq)))
    map₁, map₂ = [Vector{Int64}() for _ in 1:length(nspin), _ in 1:length(norbi)], [Vector{Int64}() for _ in 1:length(unitcell)]
    channels = Vector{Tuple{Vector{Int64},Vector{Int64}}}()
    for i in eachindex(nspin), j in eachindex(norbi)
        for k in eachindex(seq)
            if nspin[i] == seq[k][1]&&norbi[j]==seq[k][3]
                push!(map₁[i, j], k)
            end
        end
    end
    for i in 1:length(unitcell)
        for j in 1:length(cluster)
            if issubordinate(cluster.coordinates[:,j]-unitcell.coordinates[:,i], unitcell.vectors) 
                push!(map₂[i], j)
            end
        end
    end
    for i in 1:size(map₁, 1), j in 1:size(map₁, 2) 
        for u in 1:size(map₁, 1), v in 1:size(map₁, 2) 
            push!(channels, ([i, j], [u, v]))
        end
    end
    return Perioder(map₁, map₂, channels)
end

"""
    VCA{L<:AbstractLattice, G<:GFSolver} <: Frontend

Variational Cluster Approach(VCA) method for a quantum lattice system.
"""
struct VCA{L<:AbstractLattice, G<:GFSolver} <: Frontend
    unitcell::L
    cluster::L
    origigenerator::OperatorGenerator
    refergenerator::OperatorGenerator
    solver::G
    perioder::Perioder
end

"""
    VCA(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, cons::Conservation neighbors::Neighbors, m::Int=200)

Construct the Variational Cluster Approach(VCA) method for a quantum lattice system with EDSolver
"""
function VCA(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, cons::Conservation; neighbors::Neighbors, m::Int=200)
    k = EDKind(typeof(origiterms))
    table = Table(hilbert, Metric(k, hilbert))
    origibonds = bonds(cluster, neighbors)
    referbonds = filter(bond -> isintracell(bond), origibonds)
    origigenerator, refergenerator = OperatorGenerator(origiterms, origibonds, hilbert; table = table), OperatorGenerator(referterms, referbonds, hilbert; table = table)
    edsolver = EDSolver(k, refergenerator, cons, table; m = m)
    perioder = Perioder(unitcell, cluster, table)
    return VCA{typeof(unitcell), typeof(edsolver)}(unitcell, cluster, origigenerator, refergenerator, edsolver, perioder)
end
function VCA(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, cons::Conservation, rparam::Parameters; neighbors::Neighbors, m::Int=200)
    k = EDKind(typeof(origiterms))
    table = Table(hilbert, Metric(k, hilbert))
    origibonds = bonds(cluster, neighbors)
    referbonds = filter(bond -> isintracell(bond), origibonds)
    origigenerator = OperatorGenerator(origiterms, origibonds, hilbert; table = table)
    perioder = Perioder(unitcell, cluster, table)
    refergenerator = OperatorGenerator(referterms, referbonds, hilbert; table = table)
    update!(refergenerator; rparam...)
    edsolver = EDSolver(k, refergenerator, cons, table; m = m)
    return VCA{typeof(unitcell), typeof(edsolver)}(unitcell, cluster, origigenerator, refergenerator, edsolver, perioder)
end

"""
    VCAs{V<:VCA, S<:AbstractArray{V}}
    VCAs(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, cons::Conservation, params::AbstractArray{<:Parameters}; neighbors::Neighbors, m::Int=200)

Obtain a set of VCA
"""
struct VCAs{V<:VCA, S<:AbstractArray{V}}
    vcas::S
end
function VCAs(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, cons::Conservation, params::AbstractArray{<:Parameters}; neighbors::Neighbors, m::Int=200)
    k = EDKind(typeof(origiterms))
    table = Table(hilbert, Metric(k, hilbert))
    origibonds = bonds(cluster, neighbors)
    referbonds = filter(bond -> isintracell(bond), origibonds)
    origigenerator = OperatorGenerator(origiterms, origibonds, hilbert; table = table)
    perioder = Perioder(unitcell, cluster, table)
    vcas = Vector{VCA}(undef, length(params))
    for i in eachindex(params)
        refergenerator = OperatorGenerator(referterms, referbonds, hilbert; table = table)
        update!(refergenerator; params[i]...)
        edsolver = EDSolver(k, refergenerator, cons, table; m = m)
        vcas[i] = VCA{typeof(unitcell), typeof(edsolver)}(unitcell, cluster, origigenerator, refergenerator, edsolver, perioder)
    end 
    return vcas
end

"""
    origiQuadraticTerms!(om::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector,k::AbstractVector)

Calculate the Hamiltonian's quadratic terms of the original system
"""
function origiQuadraticTerms!(om::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector,k::AbstractVector)
    for i in eachindex(oops)
        @views seq₁, seq₂ = oopsseqs[i][1], oopsseqs[i][2]
        phase = exp(im*dot(k, icoordinate(oops[i])))
        om[seq₁, seq₂] += oops[i].value*phase
    end
    return om
end

"""
    CPTcore(cgfm::AbstractMatrix, vm::AbstractMatrix)

The calculating core of the CPT method
"""
function CPTcore(cgfm::AbstractMatrix, vm::AbstractMatrix)
    result = Matrix{ComplexF64}(I, size(vm)...)
    return cgfm*inv(mul!(result, vm, cgfm, -1, 1))
end

"""
    periodization(gfml::AbstractMatrix, map₂::AbstractVector, coordinates::AbstractMatrix, k::AbstractVector)

perform the periodic procedure
"""
function periodization(gfml::AbstractMatrix, map₂::AbstractVector, coordinates::AbstractMatrix, k::AbstractVector)
    N, L = length(map₂), size(coordinates,2)
    lgfm, pgfm =Matrix{ComplexF64}(undef, L, L), Matrix{ComplexF64}(undef, N, N)
    for i in 1:L, j in 1:L
        @views ra, rb = coordinates[:, i], coordinates[:, j]
        lgfm[i, j] = gfml[i,j]*exp(-im*dot(k, (ra - rb)))
    end
    for m in 1:N, n in 1:N
        @views cmap₂, hmap₂ = map₂[m], map₂[n]
        pgfm[m, n] = (N/L)*sum(lgfm[cmap₂, hmap₂])
    end
    return pgfm
end

"""
    perGreenFunction(GFm::AbstractMatrix, k::AbstractVector, perioder::Perioder, cluster::AbstractLattice)

The periodization Green function with respect to whole k-space
"""
function perGreenFunction(GFm::AbstractMatrix, k::AbstractVector, perioder::Perioder, cluster::AbstractLattice)
    gfv = Vector{Matrix{ComplexF64}}(undef, length(perioder.channels))
    gfm = Matrix{ComplexF64}(undef, length(perioder.map₁)*length(perioder.map₂), length(perioder.map₁)*length(perioder.map₂))
    for i in eachindex(perioder.channels)
        gfml = GFm[perioder.map₁[perioder.channels[i][1]...], perioder.map₁[perioder.channels[i][2]...]]
        gfv[i] = periodization(gfml, perioder.map₂, cluster.coordinates, k)
    end
    gfmm = reshape(gfv, (length(perioder.map₁), length(perioder.map₁)))
    for i in 1:length(perioder.map₁), j in 1:length(perioder.map₁)
        for u in 1:length(perioder.map₂), v in 1:length(perioder.map₂)
            gfm[(i-1)*length(perioder.map₂) + u, (j-1)*length(perioder.map₂) + v] = gfmm[i, j][u, v]
        end
    end
    return gfm
end

"""
    GreenFunctionPath(om::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector, rm::AbstractMatrix, perioder::Perioder, cluster::AbstractLattice, k_path::ReciprocalSpace, CGFm::AbstractMatrix)

Give the Green function of a certain path or area in k-space with respect to a certain energy
"""
function GreenFunctionPath(om::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector, rm::AbstractMatrix, perioder::Perioder, cluster::AbstractLattice, k_path::ReciprocalSpace, CGFm::AbstractMatrix)
    gfpath = Vector{Matrix{ComplexF64}}(undef, length(k_path))
    for i in eachindex(k_path)
        dest = copy(om)
        Vm = origiQuadraticTerms!(dest, oops, oopsseqs, k_path[i]) - rm
        GFm = CPTcore(CGFm, Vm)
        gfpath[i] = perGreenFunction(GFm, k_path[i], perioder, cluster)
    end
    return gfpath
end 

"""
    referQuadraticTerms(rops::AbstractVector, table::Table)
    
Calculate  the Hamiltonian's quadratic terms of the reference system
"""
function referQuadraticTerms(rops::AbstractVector, table::Table)
    rm = zeros(ComplexF64, length(table), length(table))
    for rop in rops 
        seq₁, seq₂ = table[rop[1].index'], table[rop[2].index]
        rm[seq₁, seq₂] += rop.value
    end
    return rm
end

function origiQuadraticTerms(oops::AbstractVector, table::Table, k::AbstractVector)
    om = zeros(ComplexF64, length(table), length(table))
    for oop in oops 
        seq₁, seq₂ = table[oop[1].index'], table[oop[2].index]
        phase = exp(im*dot(k, icoordinate(oop)))
        om[seq₁, seq₂] += oop.value*phase
    end
    return om
end

"""
    seqs(oops::AbstractVector, table::Table)

The index's sequences of original system operators
"""
function seqs(oops::AbstractVector, table::Table)
    seqs = [zeros(Int, 2) for _ in 1:length(oops)]
    for i in eachindex(oops) 
        seq₁, seq₂ = table[oops[i][1].index'], table[oops[i][2].index]
        seqs[i] = [seq₁, seq₂]
    end
    return seqs
end

"""
    singleParticleGreenFunction(vca::VCA, k_path::ReciprocalSpace, ω_range::AbstractRange)

The single particle Green function in k-ω space
"""
function singleParticleGreenFunction(vca::VCA, k_path::ReciprocalSpace, ω_range::AbstractRange)
    oops, rops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator))), filter(op -> length(op) == 2, collect(expand(vca.refergenerator)))
    oopsseqs = seqs(oops, vca.origigenerator.table)
    cgfv = clusterGreenFunctionLoop(vca.solver, ω_range)
    rm = referQuadraticTerms(rops, vca.refergenerator.table)
    om = zeros(ComplexF64, length(vca.refergenerator.table), length(vca.refergenerator.table))
    gfpv = [[Matrix{ComplexF64}(undef, length(vca.perioder.map₁)*length(vca.perioder.map₂), length(vca.perioder.map₁)*length(vca.perioder.map₂)) for _ in 1:length(k_path)] for _ in 1:length(ω_range)]
    for i in eachindex(ω_range)
        gfpv[i] = GreenFunctionPath(om, oops, oopsseqs, rm, vca.perioder, vca.cluster, k_path, cgfv[i])
    end
    return gfpv
end

"""
    Optimal{T<:VCA, P<:Parameters}

the saddle point in the reference sace
"""
struct Optimal{T<:VCA, P<:Parameters}
    optvca::T
    optparams::P
    function Optimal(optvca::VCA, optparams::Parameters)
        new{typeof(optvca),typeof(optparams)}(optvca, optparams)
    end
end

"""
    optimals(vcas::AbstractArray{<:VCA}, gps::AbstractArray{<:Real},varparams::AbstractArray{<:Parameters})

select saddle points from reference spaces
"""
function optimals(vcas::AbstractArray{<:VCA}, gps::AbstractArray{<:Real},varparams::AbstractArray{<:Parameters})
    if ndims(varparams)==2
        opts = Vector{Optimal}(undef, size(varparams, 1))
        for i in 1:size(varparams,1)
            Δgps = gps[i,:] .- maximum(gps[i,:])
            index = argmin(Δgps)
            opts[i] = Optimal(vcas[i, index], varparams[i, index])
        end
        return opts
    end
end

"""
    OrderParamters
"""
function OPintegrand(vca::VCA, bz::ReciprocalSpace, iω::Complex, sm::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector, rm::AbstractMatrix, μ::Real)
    om = zeros(ComplexF64, length(vca.refergenerator.table), length(vca.refergenerator.table))
    cgfm = clusterGreenFunction(vca.solver, iω+μ)
    intra = 0.0
    for k in bz
        vm = origiQuadraticTerms!(copy(om), oops, oopsseqs, k) - rm
        gfm = CPTcore(cgfm, vm)
        intra += (tr(sm*gfm) - tr(sm)/(iω-1.0)).re
    end
    return intra
end

function OrderParameters(opt::Optimal, hilbert::Hilbert, bz::ReciprocalSpace, term::Term, μ::Real)
    vca = opt.optvca
    term.value = convert(typeof(term.value), 1.0)
    sm = referQuadraticTerms(collect(expand(term, filter(bond -> isintracell(bond), bonds(vca.cluster, term.bondkind)), hilbert)), vca.refergenerator.table)
    oops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator)))
    oopsseqs = seqs(oops, vca.origigenerator.table)
    rm = referQuadraticTerms(filter(op -> length(op) == 2, collect(expand(vca.refergenerator))), vca.refergenerator.table)
    return abs((1/length(bz))*quadgk(x -> OPintegrand(vca, bz, x*im, sm, oops, oopsseqs, rm, μ), 0, Inf)[1]/π/length(vca.unitcell)/length(vca.cluster))
end

function clusterspectrum(cgfpathv::AbstractVector)
    A = zeros(Float64, length(cgfpathv))
    for j in eachindex(cgfpathv)
        A[j] = (-1/π) * (tr(cgfpathv[j])).im
    end
    return A
end

"""
    spectrum(gfpathv::AbstractVector)

The spectrum of the single particle Green function
"""
function spectrum(gfpathv::AbstractVector)
    A = zeros(Float64, length(gfpathv), length(gfpathv[1]))
        for i in eachindex(gfpathv)
            for j in eachindex(gfpathv[i])
                A[i, j] = (-1/π) * (tr(gfpathv[i][j])).im
            end
        end
        return A
end

"""
    GPcore(temp::AbstractMatrix, cgfm::AbstractMatrix, vm::AbstractMatrix)

The calculating core of the integrand function of the grand potential
"""
function GPcore(temp::AbstractMatrix, cgfm::AbstractMatrix, vm::AbstractMatrix)
    result = copy(temp)
    mul!(result, vm, cgfm, -1, 1)
    return log(abs(det(result)))
end

"""
    GPintegrand(solver::GFSolver, temp::AbstractMatrix, vmvec::AbstractVector, bz::AbstractVector, ω::Complex)
The integrand function of the grand potential
"""
function GPintegrand(solver::GFSolver, temp::AbstractMatrix, vmvec::AbstractVector, ω::Complex)
    intra = 0.0
    cgfm = clusterGreenFunction(solver, ω)
    for i in eachindex(vmvec)          
        intra += GPcore(temp, cgfm, vmvec[i])
    end
    return intra
end 

"""
    GrandPotential(vca::VCA, bz::AbstractVector, μ::Real) 

The grand potential for a certain VCA
"""
function GrandPotential(vca::VCA, bz::AbstractVector, μ::Real) 
    oops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator)))
    oopsseqs = seqs(oops, vca.origigenerator.table)
    rm = referQuadraticTerms(filter(op -> length(op) == 2, collect(expand(vca.refergenerator))), vca.refergenerator.table)
    om = zeros(ComplexF64, length(vca.origigenerator.table), length(vca.origigenerator.table))
    temp = Matrix{ComplexF64}(I, size(om)...)
    vmvec = Vector{Matrix{ComplexF64}}(undef,length(bz))
    for i in eachindex(bz)
        dest = copy(om)
        vmvec[i] = origiQuadraticTerms!(dest, oops, oopsseqs, bz[i]) - rm
    end
    trvm = sum([tr(vmv) for vmv in vmvec])
    gp = (vca.solver.sysvals.gsenergy + (1/length(bz))*(- quadgk(x -> GPintegrand(vca.solver, temp, vmvec, x*im+μ), 0, Inf)[1]/π + trvm.re/2))/length(vca.cluster)/length(vca.unitcell)
    return gp
end

function showgpint(vca::VCA, bz::AbstractVector, μ::Real, ω::Real)
    oops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator)))
    rm = referQuadraticTerms(filter(op -> length(op) == 2, collect(expand(vca.refergenerator))), vca.refergenerator.table)
    om = zeros(ComplexF64, length(vca.refergenerator.table), length(vca.refergenerator.table))
    temp = Matrix{ComplexF64}(I, size(om)...)
    vmvec = Vector{Matrix{ComplexF64}}(undef,length(bz))
    for i in eachindex(bz)
        vmvec[i] = origiQuadraticTerms(oops,vca.origigenerator.table, bz[i]) - rm
    end
    return GPintegrand(vca.solver, temp, vmvec, ω*im+μ)
end

function showvm(vca::VCA, bz::AbstractVector)
    oops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator)))
    oopsseqs = seqs(oops, vca.origigenerator.table)
    rm = referQuadraticTerms(filter(op -> length(op) == 2, collect(expand(vca.refergenerator))), vca.refergenerator.table)
    om = zeros(ComplexF64, length(vca.refergenerator.table), length(vca.refergenerator.table))
    vmvec = Vector{Matrix{ComplexF64}}(undef,length(bz))
    for i in eachindex(bz)
        dest = copy(om)
        vmvec[i] = origiQuadraticTerms!(dest, oops, oopsseqs, bz[i]) - rm
    end
    return vmvec
end
"""
    antiferro(wavevector::AbstractVector)

The phase of a weiss field with antiferromagnetic order
"""
function antiferro(wavevector::AbstractVector)
    function bondfunction(bond::Bond)
        return exp(im*dot(wavevector,rcoordinate(bond)))
    end
    return bondfunction
end

"""
    spiral order 
"""
function belongs(ids::AbstractVector)
    function sublattices(bond::Bond)
        if bond.points[1].site in ids
            return 1.0
        else
            return 0.0
        end
    end
    return sublattices
end
"""
    moirephase(ϕ::Real, spin::Rational)

The phase of hopping terms in moire Hubbard model
"""
function moirephase(ϕ::Real, spin::Rational)
    function spinup(bond::Bond)
        θ = azimuth(rcoordinate(bond))
        @assert any(≈(θ),(0, 2π/3, 4π/3, π/3, π, 5π/3, 2π)) "Triangle error:wrong input bond."
        any(≈(θ),(0, 2π/3, 4π/3, 2π)) && return exp(-im*ϕ)
        any(≈(θ),(π/3, π, 5π/3)) && return exp(im*ϕ)
    end
    function spindown(bond::Bond)
        θ = azimuth(rcoordinate(bond))
        @assert any(≈(θ),(0, 2π/3, 4π/3, π/3, π, 5π/3, 2π)) "Triangle error:wrong input bond."
        any(≈(θ),(0, 2π/3, 4π/3, 2π)) && return exp(im*ϕ)
        any(≈(θ),(π/3, π, 5π/3)) && return exp(-im*ϕ)
    end
    if spin == 1//2
        return spinup
    elseif spin == -1//2
        return spindown
    else
        return "error"
    end
end

function spinrotation(θ::Real, ϕ::Real)
    sin(θ)*cos(ϕ)*MatrixCoupling(:, FID, :, σ"x", :) + sin(θ)*sin(ϕ)*MatrixCoupling(:, FID, :, σ"y", :) + cos(θ)*MatrixCoupling(:, FID, :, σ"z", :)
end

"""
    spawn(numworkers::Int)

Enable multiprocess
"""
function spawn(numworkers::Int)
    np = length(workers())
    if np == 1
        addworkers = numworkers
    elseif np < numworkers
        addworkers = numworkers - np
    else
        addworkers = 0
    end
    addprocs(addworkers)
    @everywhere collect(2:numworkers+1) include("./src/parallelusing.jl")
end

"""
    saveData(data, filename::String) -> .jls

save data(e.g. a VCA data) as a jls file
"""
function saveData(data, filename::String)
    open(filename, "w") do io
        serialize(io, data)
    end
end

"""
    loadData(filename::String)

load data from a jls file
"""
function loadData(filename::String)
    data = nothing
    open(filename, "r") do io
        data = deserialize(io)
    end
    return data
end

end # module CPTVCA
