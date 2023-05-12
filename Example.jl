using Arpack: eigs
using LinearAlgebra: Eigen, dot
using SparseArrays: SparseMatrixCSC, spzeros
using QuantumLattices
using ExactDiagonalization
using BlockArrays
using Revise
using CPTVCA

lattice = Lattice([0, 0], [1, 0], [0,1], [1,1])

cluster = Cluster(lattice)

vectors = [[2,0],[0,2]]

clusters = Clusters(cluster,(1,1),vectors)


hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))
braket = (BinaryBases(8,4),BinaryBases(8,4))
table = Table(hilbert)

ops = interhopping(clusters,cluster, hilbert, Neighbors(1=>1.0),[pi/2,pi/2])

#= ops_dict = Dict{Tuple{Int,Int}, Vector{Operator}}()
hopping = Dict{Tuple{Int,Int}, Vector{Operator}}()
for op in ops
    key = getsites(op)
    if !haskey(ops_dict, key)
        ops_dict[key] = []
    end
    push!(ops_dict[key], op)
end
for i = 1:length(lattice), j = 1:length(lattice)
    key = (i, j)
    hopping[key] = []
end
for key in keys(ops_dict)
    if haskey(hopping, key)
        hopping[key] = ops_dict[key]
    end
end =#
