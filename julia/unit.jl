require("createParamFile.jl")
require("functions.jl")


println("given parameters")

g = 9.81
println("\ng = ",g)

rho = 1
println("\nrho = ",rho)

delRho = 0.5
println("\ndelta rho = ",delRho)

diameter = 0.01
println("\ndiameter = ",diameter)

sigma = 1e-4
println("\nsigma = ",sigma)

tau = 0.8
println("\ntau = ",tau)


println("\n\n##################")
println("resulting parameters")

dx, ut, c_s, dt, nu, mu, Re, Mo, Eo = getOtherParams([g,rho,delRho,diameter,sigma,tau])

println("\dx = ",dx)
println("\nterminal rise velocity = ",ut)
println("\nspeed of sound = ",c_s)
println("\ntimestep = ",dt)
println("\nnu = ",nu)
println("\nmu = ",mu)
println("\nReynolds = ",Re)
println("\nMorton = ",Mo)
println("\nEotvos = ",Eo)

createPFile(Re, Mo, Eo, tau, rho, delRho, sigma, c_s, dx, g)
