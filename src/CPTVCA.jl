module CPTVCA

using Arpack
using ExactDiagonalization
using ExactDiagonalization: matrix
using KrylovKit
using LinearAlgebra
using LinearAlgebra: lmul!
using QuantumLattices
using Serialization
using TimerOutputs

export GreenFunction, InitialState, FInitialState, SInitialState, Kryvals, LehmannGreenFunction, Sysvals, Perioder, VCA, GFSolver, EDSolver, QMCSolver
export initialstate, genkrylov, thomas, clusterGreenFunctionLoop, origiQuadraticTerms!, referQuadraticTerms, CPTcore, perGreenFunction, GreenFunctionPath, periodization, singleParticleGreenFunction, spectrum, energysurface, select, statedensity, saveData, loadData

const vcatimer = TimerOutput()
"""
    abstract type GreenFunction
"""
abstract type GreenFunction end

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
    genkrylov(matrix::AbstractMatrix, sstate::AbstractVector, m::Int)

Generate a krylov subspace with Lanczos iteration method
"""
function genkrylov(matrix::AbstractMatrix, sstate::AbstractVector, orth::ModifiedGramSchmidt, m::Int)
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
    Kryvals

The information obtained with the krylov subspace method that needed in calculating the cluster Green function
"""
struct Kryvals{T<:AbstractVector, S<:AbstractVector, R<:AbstractVector, N<:AbstractVector, P<:AbstractMatrix}
    avs::T
    bvs::S
    cvs::R
    norms::N
    projects::P
    function Kryvals(matrix::AbstractMatrix, initialstates::AbstractVector, m::Int=200)
        n = length(initialstates)
        avs, bvs, cvs, norms, projects = Vector{Vector}(undef,n), Vector{Vector}(undef,n), Vector{Vector}(undef,n), Vector{Float64}(undef,n), Matrix{Vector{Float64}}(undef,n,n)
        orth = KrylovKit.ModifiedGramSchmidt()
        for i in 1:n
            krybasis, avs[i], bvs[i], cvs[i] = genkrylov(matrix, initialstates[i], orth, m)
            norms[i] = √(initialstates[i]' * initialstates[i])
            for j in 1:n
                projects[i, j] = KrylovKit.project!(zeros(Float64, m), krybasis, initialstates[j])
            end
        end
        new{typeof(avs), typeof(bvs), typeof(cvs), typeof(norms), typeof(projects)}(avs, bvs, cvs, norms, projects)
    end
end
function initialstatesets(k::EDKind, gs::AbstractVector, target::TargetSpace, table::Table)
    sets₁, sets₂ = [Float64[] for _ in 1:length(table)], [Float64[] for _ in 1:length(table)]
    orderkeys = sort(collect(keys(table)), by = x -> table[x])
    for i in eachindex(orderkeys)
        ops₁, ops₂ = initialstate(k, 1, orderkeys[i]).ops, initialstate(k, 2, orderkeys[i]).ops
        initialstate₁, initialstate₂ = (matrix(ops₁, (target[2], target[1]), table)*gs)[:,1], (matrix(ops₂, (target[3], target[1]), table)*gs)[:,1]
        sets₁[i] = initialstate₁ 
        sets₂[i] = initialstate₂
    end
    return sets₁, sets₂
end

"""
    Sysvals{K<:EDKind}

The all information needed to calculate the Green Function of a finite size system 
"""
struct Sysvals{K<:EDKind, R<:Real}
    gsenergy::R
    kryvals₁::Kryvals
    kryvals₂::Kryvals
end
function Sysvals(k::EDKind, eigensystem::Eigen, ops::OperatorSum, target::TargetSpace, table::Table; m::Int=200)
    gse, gs = eigensystem.values[1], eigensystem.vectors[:,1]
    H₁, H₂ = matrix(ops, (target[2], target[2]), table), matrix(ops, (target[3], target[3]), table)
    sets₁, sets₂ = initialstatesets(k, gs, target, table)
    kryvals₁, kryvals₂ = Kryvals(H₁, sets₁, m), Kryvals(H₂, sets₂, m)
    return Sysvals{typeof(k), typeof(gse)}(gse, kryvals₁, kryvals₂)
end

"""
    The solver to calculate the cluster Green function
"""
abstract type GFSolver end

"""
    The ED solver to calculate the cluster Green function
"""
struct EDSolver <: GFSolver 
    sysvals::Sysvals
    function EDSolver(k::EDKind, refergenerator::OperatorGenerator, target::TargetSpace, table::Table; m::Int=200)
        rops = expand(refergenerator)
        Hₘ = matrix(rops, (target[1], target[1]), table)
        eigens = eigs(Hₘ; nev=1, which=:SR, tol=0.0,maxiter=300,  sigma=nothing, ritzvec=true, v0=[])
        eigensystem = Eigen(eigens[1], eigens[2])
        sysvals = Sysvals(k, eigensystem, rops, target, table; m=m)
        new(sysvals)
    end
end
"""
    The QMC solver to calculate the cluster Green function
"""
struct QMCSolver <: GFSolver end

"""
    Perioder
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
    VCA <: Frontend

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
    VCA(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert₁::Hilbert, hilbert₂::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, target::TargetSpace; neighbors::Neighbors, m::Int=200)

Construct the Variational Cluster Approach(VCA) method for a quantum lattice system with EDSolver
"""
function VCA(unitcell::AbstractLattice, cluster::AbstractLattice, hilbert::Hilbert, origiterms::Tuple{Vararg{Term}}, referterms::Tuple{Vararg{Term}}, target::TargetSpace; neighbors::Neighbors, m::Int=200)
    k = EDKind(typeof(origiterms))
    table = Table(hilbert, Metric(k, hilbert))
    origibonds = bonds(cluster, neighbors)
    referbonds = filter(bond -> isintracell(bond), origibonds)
    origigenerator, refergenerator = OperatorGenerator(origiterms, origibonds, hilbert; table = table), OperatorGenerator(referterms, referbonds, hilbert; table = table)
    edsolver = EDSolver(k, refergenerator, target, table; m = m)
    perioder = Perioder(unitcell, cluster, table)
    return VCA{typeof(unitcell), typeof(edsolver)}(unitcell, cluster, origigenerator, refergenerator, edsolver, perioder)
end

"""
    VCA

Construct the Variational Cluster Approach(VCA) method for a quantum lattice system with QMCSolver
"""
function VCA(qmcsolver::QMCSolver)
    #to be update
end

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

function solvercore(a::Complex, r::Complex, avs₁::AbstractVector, bvs₁::AbstractVector, cvs₁::AbstractVector, norm₁::AbstractVector, proj₁::AbstractMatrix, avs₂::AbstractVector, bvs₂::AbstractVector, cvs₂::AbstractVector, norm₂::AbstractVector, proj₂::AbstractMatrix, d::AbstractVector, m::Int, gfm::AbstractMatrix)
    sgfm = copy(gfm)
    for i in eachindex(norm₁)
        @views b₁, b₂ = a .+ bvs₁[i], r .- bvs₂[i]
        @views tmp₁, tmp₂, n₁, n₂ = thomas(avs₁[i],b₁,cvs₁[i],d,m), thomas(-avs₂[i],b₂,-cvs₂[i],d,m), norm₁[i], norm₂[i]
        for j in eachindex(norm₁)
            @views proj₁_slice, proj₂_slice = proj₁[:, j][i], proj₂[:, j][i]
            sgfm[i, j] = dot(proj₁_slice, tmp₁) * n₁ + dot(proj₂_slice, tmp₂) * n₂
        end
    end
    return sgfm
end

function clusterGreenFunction(solver::GFSolver, ω::Real, μ::Real, η::Real)
    kry₁, kry₂ = solver.sysvals.kryvals₁, solver.sysvals.kryvals₂
    avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂ = kry₁.avs, kry₁.bvs, kry₁.cvs, kry₁.norms, kry₁.projects, kry₂.avs, kry₂.bvs, kry₂.cvs, kry₂.norms, kry₂.projects
    gfm = zeros(ComplexF64, length(norm₁), length(norm₁))
    m = length(avs₁[1])
    d = [1.0;zeros(m-1)]
    a, r = ω + η*im + μ - solver.sysvals.gsenergy, ω + η*im + μ + solver.sysvals.gsenergy
    return solvercore(a, r, avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂, d, m, gfm)
end

function clusterGreenFunctionLoop(solver::GFSolver, ω_range::AbstractRange, μ::Real, η::Real)
    kry₁, kry₂ = solver.sysvals.kryvals₁, solver.sysvals.kryvals₂
    avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂ = kry₁.avs, kry₁.bvs, kry₁.cvs, kry₁.norms, kry₁.projects, kry₂.avs, kry₂.bvs, kry₂.cvs, kry₂.norms, kry₂.projects
    cgfv = [Matrix{ComplexF64}(undef, length(norm₁), length(norm₁)) for _ in 1:length(ω_range)]
    gfm = zeros(ComplexF64, length(norm₁), length(norm₁))
    m = length(avs₁[1])
    d = [1.0;zeros(m-1)]
    for i in eachindex(ω_range)
        a, r = ω_range[i] + η*im + μ - solver.sysvals.gsenergy, ω_range[i] + η*im + μ + solver.sysvals.gsenergy
        cgfv[i] = solvercore(a, r, avs₁, bvs₁, cvs₁, norm₁, proj₁, avs₂, bvs₂, cvs₂, norm₂, proj₂, d, m, gfm)
    end
    return cgfv
end
"""
    differQuadraticTerms(ogen::OperatorGenerator, rgen::OperatorGenerator, table::Table, k::AbstractVector) -> Matrix

Calculate the difference between the Hamiltonian's quadratic terms of the original system and a reference system
"""
@inline function origiQuadraticTerms!(om::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector,k::AbstractVector)
    for i in eachindex(oops)
        @views seq₁, seq₂ = oopsseqs[i][1], oopsseqs[i][2]
        phase = exp(im*dot(k, icoordinate(oops[i])))
        om[seq₁, seq₂] += oops[i].value*phase
    end
    return om
end

@inline function CPTcore(cgfm::AbstractMatrix, vm::AbstractMatrix)
    result = Matrix{ComplexF64}(I, size(vm)...)
    return cgfm*inv(mul!(result, vm, cgfm, -1, 1))
end

@inline function periodization(gfml::AbstractMatrix, map₂::AbstractVector, coordinates::AbstractMatrix, k::AbstractVector)
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
    VCAGreenFunction(vca::VCA, path::ReciprocalPath, CGFm::AbstractMatrix)

Calculate the Causal Green Function with VCA method in a certain path of a reciprocal space  
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
function GreenFunctionPath(om::AbstractMatrix, oops::AbstractVector, oopsseqs::AbstractVector, rm::AbstractMatrix, perioder::Perioder, cluster::AbstractLattice, k_path::ReciprocalPath, CGFm::AbstractMatrix)
    gfpath = Vector{Matrix{ComplexF64}}(undef, length(k_path))
    for i in eachindex(k_path)
        dest = copy(om)
        Vm = origiQuadraticTerms!(dest, oops, oopsseqs, k_path[i]) - rm
        GFm = CPTcore(CGFm, Vm)
        gfpath[i] = perGreenFunction(GFm, k_path[i], perioder, cluster)
    end
    return gfpath
end 
function referQuadraticTerms(rops::AbstractVector, table::Table)
    rm = zeros(ComplexF64, length(table), length(table))
    for rop in rops 
        seq₁, seq₂ = table[rop[1].index'], table[rop[2].index]
        rm[seq₁, seq₂] += rop.value
    end
    return rm
end

function seqs(oops::AbstractVector, table::Table)
    seqs = [zeros(Int, 2) for _ in 1:length(oops)]
    for i in eachindex(oops) 
        seq₁, seq₂ = table[oops[i][1].index'], table[oops[i][2].index]
        seqs[i] = [seq₁, seq₂]
    end
    return seqs
end

function singleParticleGreenFunction(vca::VCA, k_path::ReciprocalPath, ω_range::AbstractRange, μ::Real; η::Real = 0.05, timer=vcatimer)
    @timeit timer "SPGF" begin
    oops, rops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator))), filter(op -> length(op) == 2, collect(expand(vca.refergenerator)))
    oopsseqs = seqs(oops, vca.origigenerator.table)
    cgfv = clusterGreenFunctionLoop(vca.solver, ω_range, μ, η)
    rm = referQuadraticTerms(rops, vca.refergenerator.table)
    om = zeros(ComplexF64, length(vca.refergenerator.table), length(vca.refergenerator.table))
    gfpv = [[Matrix{ComplexF64}(undef, length(vca.perioder.map₁)*length(vca.perioder.map₂), length(vca.perioder.map₁)*length(vca.perioder.map₂)) for _ in 1:length(k_path)] for _ in 1:length(ω_range)]
    for i in eachindex(ω_range)
        gfpv[i] = GreenFunctionPath(om, oops, oopsseqs, rm, vca.perioder, vca.cluster, k_path, cgfv[i])
    end
end
    return gfpv
end

"""
    singleparticlespectrum(vca::VCA, k_path::ReciprocalPath, ω_range::AbstractRange, μ::Real)

Construct the k-ω matrix to store the data of single particle spectrum.
"""

function spectrum(gfpathv::AbstractVector)
    A = zeros(Float64, length(gfpathv), length(gfpathv[1]))
        for i in eachindex(gfpathv)
            for j in eachindex(gfpathv[i])
                A[i, j] = (-1/π) * imag(tr(gfpathv[i][j]))
            end
        end
        return A
end
function statedensity(vca::VCA, k_path::ReciprocalPath, ω_range::AbstractRange, μ::Real; η::Real=0.05)
    A = singleparticlespectrum(vca, k_path, ω_range, μ; η = η)
    S = sum(A, dims=2)
    return S
end

function energysurface(vca::VCA, reciprocalspace::ReciprocalZone, ω_range::AbstractRange, μ::Real, a::Real; η::Real=0.05)
    oops, rops = filter(op -> length(op) == 2, collect(expand(vca.origigenerator))), filter(op -> length(op) == 2, collect(expand(vca.refergenerator)))
    oopsseqs = seqs(oops, vca.origigenerator.table)
    rm = referQuadraticTerms(rops, vca.refergenerator.table)
    om = zeros(ComplexF64, length(vca.refergenerator.table), length(vca.refergenerator.table))
    es = fill(NaN64, length(reciprocalspace))
    for ω in ω_range
        cgf = clusterGreenFunction(vca.solver, ω, μ, η)
        for j in eachindex(reciprocalspace)
            if isnan(es[j])
            dest = copy(om)
            Vm = origiQuadraticTerms!(dest, oops, oopsseqs, reciprocalspace[j]) - rm
            GFm = CPTcore(cgf, Vm)
                if (-1/π)*imag(tr(perGreenFunction(GFm, reciprocalspace[j], vca.perioder, vca.cluster))) > a
                    es[j] = ω
                end
            end
        end
    end
    for e in es
        if isnan(e)
            if first(ω_range) > last(ω_range)
                e = last(ω_range) - 1 
            else
                e = last(ω_range) + 1 
            end
        end
    end
    m = reshape(es, (Int(sqrt(length(reciprocalspace))), Int(sqrt(length(reciprocalspace)))))
    reverse!(m, dims=1)
    return m
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
