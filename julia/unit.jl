using Optim

require("createParamFile.jl")
require("functions.jl")


println("given parameters")

g = 9.81
println("\ng = ",g)

diameter = 0.01
println("\ndiameter = ",diameter)

sigma = 1e-4
println("\nsigma = ",sigma)


res = optimize(targetF, [g, diameter, sigma])
show(res)


println("\n\n##################")
println("resulting parameters")

Re, Mo, Eo, dx, ut, c_s, dt, nu, mu, rho, delRho, tau = getOtherParams(g,diameter,sigma)

zwischen = targetF([g,diameter,sigma])

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
