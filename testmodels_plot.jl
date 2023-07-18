using CPTVCA
using ExactDiagonalization
using QuantumLattices

#square lattice L = (2,2)
vca = loadData("square(2,2)U4_vca.jls")
path = @rectangle_str("Γ-X-M-Γ")
k_path = selectpath(BrillouinZone(reciprocals(unitcell.vectors), (200, 200)), path)[1]
ω_range = range(-6, 6, length=600)
A = singleparticlespectrum(vca, k_path, ω_range, U.value/2)
specfig = heatmap(1:size(A)[2], ω_range, A, xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function(L=2, U=4, n=1/2)",clims=(0, 3))
N = size(A)[2]
ratios = [1, 1, sqrt(2)]
split_points = round.(Int, cumsum(ratios) ./ sum(ratios) * N)
special_points = pushfirst!(split_points,1)
tick_labels = ["Γ", "X", "M", "Γ"]
xticks!(special_points, tick_labels)
savefig("square(2,2)U4.pdf")