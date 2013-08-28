function createPFile(Re, Mo, Eo, tau, rho, delRho, sigma, c_s, dx, g)
	omega = 1/tau
	gamma = rho / (rho - delRho)

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
end

println("createParamFile wurde geladen")