using Arpack: eigs
using LinearAlgebra: Eigen
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

interlinks(cluster,clusters[8],Neighbors(1=>1.0))
filter(x->x!=cluster,clusters)[8]
typeof(interlinks(cluster,clusters[8],Neighbors(1=>1.0)))