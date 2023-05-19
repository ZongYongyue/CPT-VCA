using Revise
using QuantumLattices
using ExactDiagonalization
using CPTVCA
using Plots

lattice = Lattice([0, 0], [1, 0], [0,1], [1,1])

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

bases = BinaryBases(8, 4)

terms = (Hopping(:t, -1.0, 1), Hubbard(:U, 8.0))

vectors = [[1.0, 0.0], [0.0, 1.0]]

neighbors = Neighbors(1=>1.0)

kx_range = range(-π, π, length = 100)

ω_range = range(-1, 1, length = 100)

A_values = zeros(length(kx_range), length(ω_range))

for (i, k) in enumerate(kx_range)
    for (j, ω) in enumerate(ω_range)
        A_values[i, j] = (-1 / π) * imag(CPTGF(lattice, hilbert, terms, bases, neighbors, [k ,0], ω, vectors))
    end
end

#heatmap(kx_range, ω_range, A_values, xlabel="kₓ", ylabel="ω", color=:plasma, title="Spectral Function")

A_values

