using Optim

require("createParamFile.jl")
require("functions.jl")


println("given parameters")

g = 10.0
diameter = 1.0
sigma = 10000.0
# delRho = 0.5


res = optimize(targetF, [g, diameter, sigma], method = :nelder_mead, iterations = 10000)
show(res)

x = transform(res.minimum)

g = 		x[1]
diameter = 	x[2]
#sigma = 	1e-4 * x[3]

Re, Mo, Eo, dx, ut, c_s, dt, nu, mu, rho, delRho, tau = getOtherParams(g,diameter, sigma)

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
println("\nReynolds = ",Re)
println("\nMorton = ",Mo)
println("\nEotvos = ",Eo)

createPFile(Re, Mo, Eo, tau, rho, delRho, sigma, c_s, dx, g)
