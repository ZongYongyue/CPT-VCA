include("GPweiss.jl")
using Distributed
using Plots
spawn(11)

@time vcas = pmap(param -> VCA(unitcell, cluster, hilbert, origiterms, referterms, cons, param; neighbors=neighbors, m=200), varparams)
#saveData(vcas, "./temp/squareL4_af.jls")
#vcas = loadData("./temp/squareL4_af.jls")
#@time gps = pmap(vca -> GrandPotential(vca, rz, real(Parameters(vca.refergenerator)[:U]/2)), vcas)
@time gps = pmap(vca -> GrandPotential(vca, rz, 0), vcas)
#saveData(gps, "./temp/squareL4af_gp.jls")
#gps = loadData("./temp/squareL4af_gp.jls")

gses = [vca.solver.sysvals.gsenergy for vca in vcas]




#opts = optimals(vcas, gps, varparams)
#@time ops = pmap(opt -> OrderParameters(opt, hilbert, rz, af, opt.optparams[:U]/2), opts)