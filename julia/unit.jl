using Optim

require("createParamFile.jl")
require("functions.jl")


println("given parameters")

g = 10.0
sigma = 1e-4
diameter = 0.1

res = optimize(targetF, [g, sigma, diameter], method = :nelder_mead, iterations = 10000)
show(res)

x = transform(res.minimum)

g = 		x[1]
sigma = 	x[2]
diameter = 	x[3]

Re, Mo, Eo, dx, ut, c_s, dt, nu, mu, rho, delRho, tau = getOtherParams(g, sigma, diameter)

println("\n\n##################")
println("resulting parameters")

println("\dx = ",dx)
println("\nterminal rise velocity = ",ut)
println("\nspeed of sound = ",c_s)
println("\ntimestep = ",dt)
println("\nnu = ",nu)
println("\nmu = ",mu)
println("\nrho = ",rho)
println("\ndelta rho = ",delRho)
println("\ntau = ",tau)

println("\ng = ",g)
println("\ndiameter = ",diameter)
println("\nsigma = ",sigma)

println("\nReynolds = ",Re)
println("\nMorton = ",Mo)
println("\nEotvos = ",Eo)

createPFile(Re, Mo, Eo, tau, rho, delRho, sigma, c_s, dx, g)
