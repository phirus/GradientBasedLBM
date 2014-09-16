using Optim

require("createParamFile.jl")
require("functions.jl")
require("constants.jl")

# start values
g = 1e-4 #10.0
sigma = 1e-4
gamma = 2

res = optimize(targetF, [g, sigma], method = :nelder_mead, iterations = 10000)
show(res)

x = transform(res.minimum)

g = 		x[1]
sigma = 	x[2]

Mo, Eo, dx, c_s, dt, nu, mu, delRho, tau = getOtherParams(g, sigma)

Mo_, Eo_, dx_, c_s_, dt_, nu_, mu_, delRho_, tau_, sigma_, g_ = getParams_alter()

xi = nu * MU_RATIO
s_2 = 1/(xi/(c_s^2 * dt) + 1/2)

println("\n\n##################")
println("resulting parameters")

println("\ndx = ",dx,"\t",dx_)
println("\nspeed of sound = ",c_s,"\t",c_s_)
println("\ntimestep = ",dt,"\t",dt_)
println("\nnu = ",nu,"\t",nu_)
println("\nmu = ",mu,"\t",mu_)
println("\nrho = ",RHO_L)
println("\ndelta rho = ",delRho,"\t",delRho_)
println("\ntau = ",tau,"\t",tau_)
println("\ns_2 = ",s_2)

println("\ng = ",g,"\t",g_)
println("\nsigma = ",sigma,"\t",sigma_)

println("\nReynolds = ",REYNOLDS_MAX_INI)
println("\nMorton = ",Mo,"\t",Mo_)
println("\nEotvos = ",Eo,"\t",Eo_)

createPFile(Mo, Eo, c_s, gamma, sigma, g)
