using CPTVCA
using ExactDiagonalization
using QuantumLattices
using Plots
using TimerOutputs
const vcatimer = TimerOutput()

vca = loadData("line(1,2)U4_vca.jls")
#give a path in the reciprocal space 
path = @line_str("Γ-X") 
k_path = selectpath(BrillouinZone(reciprocals(vca.unitcell.vectors), (200,)), path)[1]
#give the energy range
ω_range = range(-6, 6, length=400)
#caculate the spectral function of the VCA green function with half filling
A = singleparticlespectrum(vca, k_path, ω_range, Parameters(vca.origigenerator)[:U]/2, timer=vcatimer)
show(vcatimer)
#=
──────────────────────────────────────────────────────────────────────
Time                    Allocations      
───────────────────────   ────────────────────────
Tot / % measured:       10.0s /  96.1%           8.31GiB /  98.9%    

Section     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────
spectrum         1    9.62s  100.0%   9.62s   8.22GiB  100.0%  8.22GiB
CGFv           1    7.03s   73.0%   7.03s   4.05GiB   49.2%  4.05GiB
GFpathv        1    2.60s   27.0%   2.60s   4.17GiB   50.8%  4.17GiB
──────────────────────────────────────────────────────────────────────
=#
#a = select(A, (:,:))
#f = plot(k_path, ω_range, a; xlabel="k", ylabel="ω", color=:jet1, title="Spectral Function(L=2, U=4, n=1/2)",clims=(0, 3))
#savefig("tline(1,2)U4.pdf")