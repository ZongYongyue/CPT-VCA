module QuantumClusterTheory

using Arpack
using ExactDiagonalization
using KrylovKit
using LinearAlgebra
using QuantumLattices

export GreenFunction, InitialState, FInitialState, SInitialState, Kryvals, LehmannGreenFunction, Sysvals
export initialstate, clusterGreenFunction

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
        new(sign,key)
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
    Kryvals{M<:AbstractMatrix, R<:Real, V<:AbstractVector}

The information obtained with the krylov subspace method that needed in calculating the cluster Green function
"""
struct Kryvals{M<:AbstractMatrix, R<:Real, V<:AbstractVector}
    tridimatrix::M
    norm::R
    projectvector::V
    function Kryvals(matrix::AbstractMatrix, sstate::AbstractVector, m::Int=200)
        krybasis, tridimatrix = genkrylov(matrix, sstate, m)
        projectvector = KrylovKit.project!(zeros(Complex, m), krybasis, sstate)
        norm = √(sstate'*sstate)
        new{typeof(tridimatrix), typeof(norm), typeof(projectvector)}(tridimatrix, norm, projectvector)
    end
end

"""
    genkrylov(matrix::AbstractMatrix, sstate::AbstractVector, m::Int)

Generate a krylov subspace with Lanczos iteration method
"""
function genkrylov(matrix::AbstractMatrix, sstate::AbstractVector, m::Int)
    matrix = Hermitian(matrix)
    orth = KrylovKit.ModifiedGramSchmidt()
    iterator = LanczosIterator(matrix, sstate, orth)
    factorization = KrylovKit.initialize(iterator)
    for _ in 1:m-1
        KrylovKit.expand!(iterator, factorization)
    end
    basis_vectors = basis(factorization)
    T = rayleighquotient(factorization)
    return basis_vectors, T
end

"""
    LehmannGreenFunction{R<:Real, I<:Int, S<:Kryvals} <: GreenFunction

The minimum element of a Green function in Lehmann representation
"""
struct LehmannGreenFunction{R<:Real, I<:Int, S<:Kryvals} <: GreenFunction
    GSEnergy::R
    sign::I
    kryvals_l::S
    kryvals_r::S
    function LehmannGreenFunction(gse::Real, sign::Int, kryvals_l::Kryvals, kryvals_r::Kryvals)
        new{typeof(gse), typeof(sign), typeof(kryvals_l)}(gse, sign, kryvals_l, kryvals_r)
    end
end
function (gf::LehmannGreenFunction)(ω::Complex, μ::Real)
    Im = Matrix{Complex}(I, size(gf.kryvals_l.tridimatrix,2), size(gf.kryvals_l.tridimatrix,2))
    if gf.sign == 1
        lgf = dot(gf.kryvals_l.projectvector, inv((ω + μ - gf.GSEnergy)*Im + Matrix(gf.kryvals_r.tridimatrix))[:, 1])*gf.kryvals_r.norm
    elseif gf.sign == 2
        lgf = dot(gf.kryvals_l.projectvector, inv((ω + μ + gf.GSEnergy)*Im - Matrix(gf.kryvals_r.tridimatrix))[:, 1])*gf.kryvals_r.norm
    end
    return lgf
end

"""
    Sysvals{K<:EDKind, R<:Real, S<:Kryvals}

The all information needed to calculate the Green Function of a finite size system 
"""
struct Sysvals{K<:EDKind, R<:Real, S<:Kryvals}
    GSEnergy::R
    setkryvals₁::Vector{S}
    setkryvals₂::Vector{S}
end
function Sysvals(k::EDKind, gs::AbstractVector, gse::Real, gen::OperatorGenerator, target::TargetSpace, hilbert::Hilbert; m::Int=200)
    table = Table(hilbert, Metric(k, hilbert))
    setkryvals₁, setkryvals₂ = (Vector{Kryvals}(), Vector{Kryvals}())
    H₁, H₂ = matrix(expand(gen), (target[2], target[2]), table), matrix(expand(gen), (target[3], target[3]), table)
    orderkeys = sort(collect(keys(table)), by = x -> table[x])
    for key in orderkeys
        ops₁, ops₂ = initialstate(k, 1, key).ops, initialstate(k, 2, key).ops
        sstate₁, sstate₂ = (matrix(ops₁, (target[2], target[1]), table)*gs)[:,1], (matrix(ops₂, (target[3], target[1]), table)*gs)[:,1]
        push!(setkryvals₁, Kryvals(H₁, sstate₁, m))
        push!(setkryvals₂, Kryvals(H₂, sstate₂, m))
    end
    return Sysvals{typeof(k), typeof(gse), eltype(setkryvals₁)}(gse, setkryvals₁, setkryvals₂)
end

"""
    clusterGreenFunction(sys::Sysvals, ω::Complex, μ::Real)

Calculate the cluster Green function
"""
function clusterGreenFunction(sys::Sysvals, ω::Complex, μ::Real)
    cgf = zeros(Complex, length(sys.setkryvals₁), length(sys.setkryvals₁))
    for i in 1:length(sys.setkryvals₁), j in 1:length(sys.setkryvals₁)
        gf₂ = LehmannGreenFunction(sys.GSEnergy, 2, sys.setkryvals₂[i], sys.setkryvals₂[j])
        gf₁ = LehmannGreenFunction(sys.GSEnergy, 1, sys.setkryvals₁[j], sys.setkryvals₁[i])
        cgf[i, j] += gf₂(ω, μ) + gf₁(ω, μ)
    end
    return cgf
end


end # module
