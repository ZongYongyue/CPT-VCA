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

clusters = Clusters(cluster,(1,2),vectors)

interlinks(cluster,clusters[8],Neighbors(1=>1.0,3=>√(1.5^2+0.5^2)))
filter(x->x!=cluster,clusters)[8]
#typeof(interlinks(cluster,clusters[8],Neighbors(1=>1.0)))

#求集团间跃迁，希尔伯特空间是一个集团内所有格点的希尔伯特空间
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))
braket = (BinaryBases(8,4),BinaryBases(8,4))
table = Table(hilbert)
op1 = CompositeIndex(Index(1,FID{:f}(1, 1//2, 2)),[0,0],[0,0])*CompositeIndex(Index(1,FID{:f}(1, 1//2, 1)),[1,0],[1,0])
op2 = CompositeIndex(Index(1,FID{:f}(1, 1//2, 2)),[0,1],[0,1])*CompositeIndex(Index(1,FID{:f}(1, 1//2, 1)),[1,1],[1,1])
ops = op1+op1'+op2+op2'

opsm = matrix(ops,braket,table)

ed = ED(lattice,hilbert,(Hopping(:t,-1,1),),TargetSpace(braket[1]))
ed.H.operators.alterops.:t

interlinks(clusters[2],cluster,Neighbors(1=>1.0))
interlinks(clusters[4],cluster,Neighbors(1=>1.0))
interlinks(clusters[6],cluster,Neighbors(1=>1.0))
interlinks(clusters[8],cluster,Neighbors(1=>1.0))

-Rational((hilbert[1].nspin-1)/2) : Rational((hilbert[1].nspin-1)/2)
neighbors = Neighbors(1=>1.0)
clusters
ops = interhopping(clusters,cluster, hilbert,Neighbors(1=>1.0),[pi/2,pi/2])
#surroundings = filter(x->x!=cluster,clusters)
#= enumerate(surroundings)
for (index, value) in enumerate(surroundings)
    println("$index $value")
end =#
#= v = surroundings[2]
seqs = interlinks(v, cluster, neighbors)

cluster[:,seqs[1][2]] =#
#CompositeIndex(Index(seqs[1][2],FID{:f}(1, 1//2, 2)),cluster[:,seqs[1][2]], [0,0])
#CompositeIndex(Index(seqs[1][2],FID{:f}(1, 1//2, 2)), cluster[:,seqs[1][2]], [0.0,0.0])
#=
cluster[:,seqs[1][2]]
typeof(cluster[:,seqs[1][2]])
typeof([0,0])
=#
#CompositeIndex(Index(1, FID{:f}(1, 0, 2)), [0.5, 0.0], [0.0, 0.0])

matrix(ops,braket,table)