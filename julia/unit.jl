using Optim

println("##################\n")

require("createParamFile.jl")
require("functions.jl")
require("constants.jl")

println("\n##################")

Mo, Eo, dx, c_s, dt, nu, mu, delRho, tau, sigma, g = getParams_alter()

xi = nu * MU_RATIO
s_2 = 1/(xi/(c_s^2 * dt) + 1/2)


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
