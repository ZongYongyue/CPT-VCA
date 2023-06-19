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
sc_d = Pairing(:sc_d, Δ, 1, Coupling(Index(:, FID(:, 1//2, 1)), Index(:, FID(:, -1//2, 1))); amplitude = amp_sc)
terms_ref = (t_ref, μ_ref, U, AF_up, AF_down, sc_d)



lattice = Lattice([0, 0], [0.2, 0], [1, 0], [1.2, 0]; vectors=[[2, 0]])
lattice = Lattice([0, 0], [1, 0]; vectors=[[2, 0]])
lattice = Lattice([0, 0], [1, 0], [1, 0], [1, 1]; vectors=[[2, 0], [0,2]])
lattice = Lattice([0.0, 0.0], [1.0, 0.0], [1/2, √3/2]; vectors = [[3/2, √3/2], [-3/2, √3/2]])
labonds = bonds(lattice, Neighbors(0=>0, 1=>1))
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))


#生成nambu表示下的非频率依赖项的矩阵V
function V_matrix(terms::Tuple{Vararg{Term}}, terms_ref::Tuple{Vararg{Term}}, bonds, hilbert::Hilbert, k)
    table = Table(hilbert)
    om, rm = zeros(Complex, length(table), length(table)), zeros(Complex, length(table), length(table))
    refbonds = filter(bond -> isintracell(bond), bonds)
    origops, refops = expand(OperatorGenerator(terms, bonds, hilbert)), expand(OperatorGenerator(terms_ref, refbonds, hilbert))
    origops_dict = Dict(key => origops.contents[key] for key in keys(origops.contents) if length(origops.contents[key]) == 2)
    refops_dict = Dict(key => refops.contents[key] for key in keys(refops.contents) if length(refops.contents[key]) == 2)
    for (key, origop) in origops_dict
        seq₁, seq₂ = table[origop[1].index], table[origop[2].index]
        phase = isapprox(norm(icoordinate(origop)), 0.0) ? one(eltype(om)) : convert(eltype(om), 2*exp(im*dot(k, icoordinate(origop))))
        om[seq₁, seq₂] += origop.value*phase
    end
    for (key, refop) in refops_dict
        seq₁, seq₂ = table[refop[1].index], table[refop[2].index]
        rm[seq₁, seq₂] += refop.value
    end
    return om - rm
end 

V_matrix(terms, terms_ref, labonds, hilbert, [π/3, π/2])
