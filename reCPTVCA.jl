using Arpack: eigs
using LinearAlgebra
using SparseArrays: SparseMatrixCSC, spzeros, spdiagm
using QuantumLattices
using QuantumLattices: expand
using ExactDiagonalization
using ExactDiagonalization: eigen
using BlockArrays
using IterativeSolvers: minres, gmres
using StaticArrays: SVector
using KrylovKit
using CPTVCA

#original systerm terms
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 4.0)
μ = Onsite(:μ, U.value/2)

terms = (t, U, μ)

#reference systerm terms
t_ref_value = -1.1
μ_ref_value = 2.1
M = 0.1
Δ = 1
function amp_AF_up(bond::Bond)
    r = bond.points[1].rcoordinate
    return -exp(dot([π, π], r))
end
function amp_AF_down(bond::Bond)
    r = bond.points[1].rcoordinate
    return exp(dot([π, π], r))
end
function amp_sc(bond::Bond)
    θ = azimuth(rcoordinate(bond))
    any(≈(θ),(0, π)) && return 1
    any(≈(θ),(π/2, 3π/2)) && return -1
end
t_ref = Hopping(:t, t_ref_value, 1)
μ_ref = Onsite(:μ, μ_ref_value)
AF_up = Onsite(:AF_up, M, Coupling(Index(:, FID(:, 1//2, :)), Index(:, FID(:, 1//2, :))); amplitude = amp_AF_up)
AF_down = Onsite(:AF_down, M, Coupling(Index(:, FID(:, -1//2, :)), Index(:, FID(:, -1//2, :))); amplitude = amp_AF_down)
#sc_d = Pairing(:sc_d, Δ, 1, Coupling(Index(:, FID(:, 1//2, 1)), Index(:, FID(:, -1//2, 1))); amplitude = amp_sc)
terms_ref = (t_ref, μ_ref, U, AF_up, AF_down)



lattice = Lattice([0, 0], [0.2, 0], [1, 0], [1.2, 0]; vectors=[[2, 0]])
lattice = Lattice([0, 0], [1, 0]; vectors=[[2, 0]])
lattice = Lattice([0, 0], [1, 0], [1, 0], [1, 1]; vectors=[[2, 0], [0,2]])
lattice = Lattice([0.0, 0.0], [1.0, 0.0], [1/2, √3/2]; vectors = [[3/2, √3/2], [-3/2, √3/2]])
labonds = bonds(lattice, Neighbors(0=>0, 1=>1))
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))


#= 
#生成nambu表示下的非频率依赖项的矩阵V
function V_matrix(terms::Tuple{Vararg{Term}}, terms_ref::Tuple{Vararg{Term}}, bonds, hilbert::Hilbert, k)
    table = Table(hilbert)
    om, rm = zeros(Complex, length(table), length(table)), zeros(Complex, length(table), length(table))
    refbonds = filter(bond -> isintracell(bond), bonds)
    origops = filter(op -> length(op) == 2, collect(expand(OperatorGenerator(terms, bonds, hilbert))))
    refops = filter(op -> length(op) == 2,  collect(expand(OperatorGenerator(terms_ref, refbonds, hilbert))))
    for origop in origops
        seq₁, seq₂ = table[origop[1].index], table[origop[2].index]
        phase = isapprox(norm(icoordinate(origop)), 0.0) ? one(eltype(om)) : convert(eltype(om), 2*exp(im*dot(k, icoordinate(origop))))
        om[seq₁, seq₂] += origop.value*phase
    end
    for refop in refops
        seq₁, seq₂ = table[refop[1].index], table[refop[2].index]
        rm[seq₁, seq₂] += refop.value
    end
    return om - rm
end 
=#
#V_matrix(terms, terms_ref, labonds, hilbert, [π/3, π/2])

struct OHMatrixRepresentation{V, T} <: MatrixRepresentation
    k::V
    table::T
end

function (mr::OHMatrixRepresentation)(oops::OperatorSum)
    om = zeros(Complex, length(mr.table), length(mr.table))
    for oop in oops
        if length(oop) == 2
            add!(om, mr, oop)
        end
    end
    return om
end
function add!(om::AbstractMatrix, mr::OHMatrixRepresentation, oop::Operator)
    seq₁, seq₂ =  mr.table[oop[1].index], mr.table[oop[2].index]
    phase = isapprox(norm(icoordinate(oop)), 0.0) ? one(eltype(om)) : convert(eltype(om), 2*exp(im*dot(mr.k, icoordinate(oop))))
    om[seq₁, seq₂] += oop.value*phase
    return om
end

struct RHMatrixRepresentation{P, T} <: MatrixRepresentation
    parameters::P
    table::T
end
function (mr::RHMatrixRepresentation)(genrops::OperatorGenerator)
    rm = zeros(Complex, length(mr.table), length(mr.table))
    gen = update!(genrops; mr.parameters...)
    rops = filter(op -> length(op) == 2,  collect(expand(gen)))
    for rop in rops 
        if length(rop) == 2
            seq₁, seq₂ =  mr.table[rop[1].index], mr.table[rop[2].index]
            rm[seq₁, seq₂] += rop.value
        end
    end
    return rm
end

struct VCAinfo
    oops::OperatorSum
    genrops::OperatorGenerator
    table::Table
end
function VCA(info::VCAinfo, k_range::StepRangeLen)
    table = info.table
    oops = info.oops
    genrops = info.genrops
    rh = RHMatrixRepresentation(parameters, table)
    for k in k_range
        oh = OHMatrixRepresentation(k, table)
        vm = oh(oops) - rh(genrops)
    end
    return 
end

#= 
oh = OHMatrixRepresentation([π/3, π/2], Table(hilbert))
origops = expand(OperatorGenerator(terms, labonds, hilbert))
oh(origops)
rh = RHMatrixRepresentation((t = -1.2, μ = 2.1, U = 4.0, AF_up = 0.1, AF_down = 0.1, sc_d = 1), Table(hilbert))
refbonds = filter(bond -> isintracell(bond), labonds)
genrops = OperatorGenerator(terms_ref, refbonds, hilbert)
Parameters(genrops)
rh(genrops)
=#