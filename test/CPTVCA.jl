using CPTVCA
using Arpack
using ExactDiagonalization
using KrylovKit
using LinearAlgebra
using QuantumLattices
using Test

@testset "clusterGreenFunction" begin
    @test clusterGreenFunction(sys, 1+0.05im, 2) == 1
end