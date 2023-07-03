using QuantumLattices
using ExactDiagonalization
using KrylovKit
using LinearAlgebra
using Revise


abstract type GreenFunction end

abstract type OnsiteStartStates end

struct Kryvals{M<:AbstractMatrix, R<:Real, V<:AbstractVector}
    tridimatrix::M
    norm::R
    projectvector::V
    function Kryvals(matrix::AbstractMatrix, sstate::AbstractVector, m::Int=200)
        krybasis, tridimatrix = genkrylov(matrix, sstate, m)
        projectvector = KrylovKit.project!(zeros(Complex, m), krybasis, sstate)
        norm = √(sstate'*sstate)
        new{typeof(tridimatrix), typeof(norm), typeof{projectvector}}(tridimatrix, norm, projectvector)
    end
end

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


struct FOnsiteStartStates{I<:Int, O<:OperatorSum} <: OnsiteStartStates
    id::I
    contents::Vector{O}
    function FOnsiteStartStates(id::Int, site::Int, orderkeys::AbstractVector)
        contents = Vector{OperatorSum{Operator}}()
        onsitekeys = filter(key -> key[2] == site, orderkeys)
        for key in onsitekeys
            ops = Operators(Operator(1, Index(key[2], FID{:f}(key[3], key[1], id))))
            push!(contents, ops)
        end
        new{typeof(id), eltype(contents)}(id, contents)
    end
end

struct SOnsiteStartStates{I<:Int, O<:OperatorSum} <: OnsiteStartStates
    id::I
    contents::Vector{O}
end

struct LehmannGreenFunction{R<:Real, I<:Int, S<:Kryvals} 
    GSEnergy::R
    id::I
    kryvals_l::S
    kryvals_r::S
end

function (gf::LehmannGreenFunction)(ω::Complex, μ::Real)
    Im = Matrix{Complex}(I, length(gf.kryvals_l.tridimatrix), length(gf.kryvals_l.tridimatrix))
    if gf.id == 1
        lgf = dot(gf.kryvals_l.projectvector, inv((ω + μ - GSEnergy)*Im + Matrix(gf.kryvals_r.tridimatrix))[:, 1])*gf.kryvals_r.norm
    else if gf.id == 2
        lgf = dot(gf.kryvals_l.projectvector, inv((ω + μ + GSEnergy)*Im - Matrix(gf.kryvals_r.tridimatrix))[:, 1])*gf.kryvals_r.norm
    end
    return lgf
end

struct Sysvals{K<:EDKind, R<:Real, S<:Kryvals}
    GSEnergy::R
    setkryvals₁::Vector{S}
    setkryvals₂::Vector{S}
end
function Sysvals(k::EDKind{:FED}, gs::AbstractVector, gse::Real, gen::OperatorGenerator, target::TargetSpace, hilbert::Hilbert; m::Int=200) 
    table = Table(hilbert, Metric(k, hilbert))
    setkryvals₁, setkryvals₂ = fill(Vector{Kryvals}(), 2)
    H₁, H₂ = matrix(expand(gen), (target[2], target[2]), table), matrix(expand(gen), (target[3], target[3]), table)
    orderkeys = sort(collect(keys(table)), by = x -> table[x])
    for i in 1:length(hilbert)
        for ops in FOnsiteStartStates(1, i, orderkeys).contents
            sstate₁ = (matrix(ops, (target[2], target[1]), table)*gs)[:,1]
            push!(setkryvals₁, Kryvals(H₁, sstate₁, m))
        end
        for ops in FOnsiteStartStates(2, i, orderkeys).contents
            sstate₂ = (matrix(ops, (target[3], target[1]), table)*gs)[:,1]
            push!(setkryvals₂, Kryvals(H₂, sstate₂, m))
        end
    end
    return Sysvals{typeof(k), typeof{gse}, eltype{setkryvals₁}}(gse, setkryvals₁, setkryvals₂)
end

function Sysvals(k::EDKind{:SED}, gs::AbstractVector, gse::Real, gen::OperatorGenerator, target::TargetSpace, hilbert::Hilbert; m::Int=200) 
end

function ClusterGreenFunction(sys::Sysvals, ω::Complex, μ::Real)
    cgf = zeros(Complex, length(sys.setkryvals₁), length(sys.setkryvals₁))
    for i in 1:length(sys.setkryvals₁), j in 1:length(sys.setkryvals₁)
        sys.setkryvals₁[i], sys.setkryvals₁[j], sys.setkryvals₂[i], sys.setkryvals₂[j]
        gf₂ = LehmannGreenFunction(gse, 2, sys.setkryvals₂[i], sys.setkryvals₂[j])
        gf₁ = LehmannGreenFunction(gse, 1, sys.setkryvals₁[j], sys.setkryvals₁[i])
        cgf[i, j] += gf₂(ω, μ) + gf₁(ω, μ)
    end
    return cgf
end





hilbert_F = Hilbert(site=>Fock{:f}(1, 2) for site=1:2)
hilbert_S = Hilbert(site=>Spin{1//2}() for site=1:2)
st = Table(hilbert_S)
ft = Table(hilbert_F, Metric(EDKind{:FED}(),hilbert_F))

orderkeys = sort(collect(keys(ft)), by = x -> ft[x])
a = FOnsiteStartStates(1,1,orderkeys)

