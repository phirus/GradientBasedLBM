println("given parameters")

g = 9.81
println("\ng = ",g)

rho = 1
println("\nrho = ",rho)

delRho = 0.999
println("\ndelta rho = ",delRho)

diameter = 0.01
println("\ndiameter = ",diameter)

sigma = 1e-4
println("\nsigma = ",sigma)

tau = 0.8
println("\ntau = ",tau)


println("\n\n##################")
println("resulting parameters")
dx = diameter/40
println("\dx = ",dx)

ut = sqrt(2.14*sigma /(rho*diameter) + 0.505 *g * diameter)
println("\nterminal rise velocity = ",ut)

c_s = 10*ut;
println("\nspeed of sound = ",c_s)

dt = dx / (sqrt(3)*c_s)
println("\ntimestep = ",dt)

nu = c_s^2 * dt *(tau - 0.5)
println("\nnu = ",nu)
mu = rho * nu;
println("\nmu = ",mu)

Re = rho * diameter * ut / mu
Mo = g * mu^4 * delRho/(rho^2 * sigma^3)
Eo = (g*delRho*diameter^2)/sigma


println("\nReynolds = ",Re)
println("\nMorton = ",Mo)
println("\nEotvos = ",Eo)

# Pre-Output
omega = 1/tau
gamma = rho / (rho - delRho)

# Output
stream = open("parameterFile",true,true,true,true,false)
write(stream,"# Parameterset yielding")
write(stream,"\n# Reynolds = ")
write(stream, string(Re))
write(stream,"\n# Morton = ")
write(stream, string(Mo))
write(stream,"\n# Eotvos = ")
write(stream, string(Eo))

write(stream, "\n\n# relaxations parameters, no unit, default [1]")
write(stream,"\nomega_red = ")
write(stream, string(omega))
write(stream,"\nomega_blue = ")
write(stream, string(omega))

write(stream, "\n\n# density of the denser phase / kg m^-3, default [1]")
write(stream,"\nrho_red = ")
write(stream, string(rho))

write(stream, "\n\n# density ratio, default [1000]")
write(stream,"\ngamma = ")
write(stream, string(gamma))

write(stream, "\n\n# alpha, default [0.2] \nalpha_blue = 0.2")
write(stream, "\n\n# delta, default [0.1] \ndelta = 0.1 ")
write(stream, "\n\n# beta, default [0.99] \nbeta = 0.99")

write(stream, "\n\n# surface tension /?, default [1e-4]")
write(stream,"\nsigma = ")
write(stream, string(sigma))

write(stream, "\n\n# speed of sound / m s^-1, default [1484]")
write(stream,"\nc_s = ")
write(stream, string(c_s))

write(stream, "\n\n# size of a spacial step / m, default [0.001]")
write(stream,"\ndx = ")
write(stream, string(dx))


write(stream, "\n\n# gravitational constant / m s^-1 default [9.81]")
write(stream,"\ng = ")
write(stream, string(g))

write(stream,"\n")
close(stream)
