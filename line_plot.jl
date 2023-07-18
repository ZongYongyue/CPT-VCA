using CPTVCA
using ExactDiagonalization
using QuantumLattices

vca = loadData("line(1,2)U4_vca.jls")
#give a path in the reciprocal space 
path = @line_str("Γ-X") 
k_path = selectpath(BrillouinZone(reciprocals(unitcell.vectors), (400,)), path)[1]
#give the energy range
ω_range = range(-6, 6, length=600)
#caculate the spectral function of the VCA green function with half filling
A = singleparticlespectrum(vca, k_path, ω_range, U.value/2)
#plot the the spectral function
specfig = heatmap(1:size(A)[2], ω_range, A, xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function(L=2, U=4, n=1/2)",clims=(0, 3))
N = size(A)[2]
ratios =[1]
split_points = round.(Int, cumsum(ratios) ./ sum(ratios) * N)
special_points = pushfirst!(split_points,1)
tick_labels = ["Γ", "X"]
xticks!(special_points, tick_labels)
savefig("line(1,2)U4.pdf")
