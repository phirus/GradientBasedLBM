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

xi = nu * MU_RATIO
s_2 = 1/(xi/(c_s^2 * dt) + 1/2)

println("\n\n##################")
println("resulting parameters")

println("\ndx = ",dx)
println("\nspeed of sound = ",c_s)
println("\ntimestep = ",dt)
println("\nnu = ",nu)
println("\nmu = ",mu)
println("\nrho = ",RHO_L)
println("\ndelta rho = ",delRho)
println("\ntau = ",tau)
println("\ns_2 = ",s_2)

println("\ng = ",g)
println("\nsigma = ",sigma)

println("\nReynolds = ",REYNOLDS_MAX_INI)
println("\nMorton = ",Mo)
println("\nEotvos = ",Eo)

createPFile(Mo, Eo, c_s, gamma, sigma, g)
