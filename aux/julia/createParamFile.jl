function createPFile(Mo, Eo, c_s, sigma, g)

	stream = open("preprocessFile",true,true,true,true,false)
	write(stream,"# Parameterset yielding")
	write(stream,"\n Reynolds = ")
	write(stream, string(REYNOLDS_MAX_INI))
	write(stream,"\n Morton = ")
	write(stream, string(Mo))
	write(stream,"\n Eotvos = ")
	write(stream, string(Eo))

	write(stream, "\n\n# resolution / Cells per Bubble diameter, default [40]")
	write(stream,"\nresolution = ")
	write(stream, string(RESOLUTION))

	write(stream, "\n\n# density of the denser phase / kg m^-3, default [1000]")
	write(stream,"\nrho_l = ")
	write(stream, string(RHO_L))

	write(stream, "\n\n# density ratio, default [5]")
	write(stream,"\ngamma = ")
	write(stream, string(GAMMA))

	# write(stream, "\n\n# alpha, default [0.2] \nalpha_blue = 0.2")
	# write(stream, "\n\n# delta, default [0.1] \ndelta = 0.1 ")
	# write(stream, "\n\n# beta, default [0.99] \nbeta = 0.99")

	write(stream, "\n\n# speed of sound / m s^-1, default [10]")
	write(stream,"\nc_s = ")
	write(stream, string(c_s))

	write(stream, "\n\n# surface tension /?, default [1e-4]")
	write(stream,"\nsigma = ")
	write(stream, string(sigma))

	write(stream, "\n\n# gravitational constant / m s^-1, default [10]")
	write(stream,"\ng = ")
	write(stream, string(g))

	write(stream, "\n\n# bulk viscosity ratio mu'/mu  , default [2]")
	write(stream,"\nmu_ratio = ")
	write(stream, string(MU_RATIO))

	write(stream, "\n\n# s_3 , default [1]")
	write(stream,"\ns_3 = ")
	write(stream, string(S_3))

	write(stream, "\n\n# s_5 , default [1]")
	write(stream,"\ns_5 = ")
	write(stream, string(S_5))

	write(stream, "\n\n# refine factor, default [1.1]")
	write(stream,"\nfactor = ")
	write(stream, string(REFINE_FACTOR))

	write(stream, "\n\n# maximum number of time steps, default [1e5]")
	write(stream,"\nmax_steps = ")
	write(stream, string(MAX_ITER))

	write(stream, "\n\n# output Settings")
	write(stream,"\ntechplot_interval = ")
	write(stream, string(TECH_PLOT))
	write(stream,"\nrestart_interval = ")
	write(stream, string(RESTART))

	write(stream, "\n\n# termination criterion")
	write(stream,"\nresi_Re_rel = ")
	write(stream, string(RESIDUAL_RE))

	write(stream,"\n")
	close(stream)
end

println("createParamFile wurde geladen")