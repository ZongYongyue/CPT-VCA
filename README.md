# CPTVCA
Julia package for the Cluster Perturbation Theory and Variational Cluster Approach in Condensed Matter Physics.

## Introduction 
The Julia package CPTVCA, short for Cluster Perturbation Theory and Variational Cluster Approach, is a method to caculate the single particle green function of a quantum lattice systerm and the related phsical quantities based on QuantumLattices.jl and ExactDiagonalization.jl

## Main functions
interclusterbonds : obtain inter-cluster bonds between the given cluster and its surrounding clusters
interhoppingm : obtain a inter-cluster hopping matrix
CGF : obtain a cluster green function
CPTGF : obtaion a CPT green function
CPTspec : obtain the spectral function of a CPT green function 
(to be updated ... )

## Dependencies
In Julia v1.8+, please type \] in the REPL to use the package mode, then type this command:

```

pkg> add QuantumLattices
pkg> add ExactDiagonalization
pkg> add Arpack
pkg> add BlockArrays
pkg> add IterativeSolvers
pkg> add LinearAlgebra
pkg> add SparseArrays
pkg> add StaticArrays

```

## Author and Contact
Zong Yongyue (Department of Physics, Nanjing University)
zongyyphy@qq.com


