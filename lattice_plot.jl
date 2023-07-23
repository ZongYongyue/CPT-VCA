using QuantumLattices
using Plots
using LinearAlgebra
#plot lattice
#honeycomb L=6
unitcell = Lattice([0,0],[√3/2,1/2]; vectors = [[√3/2,3/2],[√3,0]])
cluster = Lattice([0,0],[√3/2,1/2],[√3,0],[√3,-1],[√3/2,-3/2],[0,-1]; vectors = [[3√3/2,3/2],[3√3/2,-3/2]])
#square L=(2, 3)
unitcell = Lattice([0, 0]; vectors=[[1, 0], [0, 1]])
cluster = Lattice([0, 0], [1, 0], [0, 1], [1, 1], [0, 2], [1, 2]; vectors = [[2, 0], [0, 3]])
neighbors = Neighbors(0=>0.0, 1=>1.0, 2=>√2)

plot(cluster,neighbors)

