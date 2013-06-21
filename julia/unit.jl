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


# Output
stream = open("parameterFile",true,true,true,true,false)
write(stream,"# Parameterset yielding")
write(stream,"\n# Reynolds = ")
write(stream, string(Re))
write(stream,"\n# Morton = ")
write(stream, string(Mo))
write(stream,"\n# Eotvos = ")
write(stream, string(Eo))

write(stream,"\nomega_red = ")
write(stream, string(omega))
write(stream,"\nomega_blue = ")
write(stream, string(omega))
write(stream,"\nrho_red = ")
write(stream, string(rho))


write(stream,"\n")
close(stream)
