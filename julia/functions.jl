function getOtherParams(g,diameter,sigma)
	rho = 		1.0
	delRho = 	0.5
	tau = 		0.8

	dx = diameter/40

	tmp = 2.14 * sigma / (rho * diameter) + 0.505 *g * diameter
	if (tmp < 0 ) 
		tmp = 0
	end

	ut = sqrt( tmp )
	c_s = 10*ut;
	dt = dx / (sqrt(3)*c_s)
	nu = c_s^2 * dt *(tau - 0.5)
	mu = rho * nu;
	Re = rho * diameter * ut / mu
	Mo = g * mu^4 * delRho/(rho^2 * sigma^3)
	Eo = (g*delRho*diameter^2)/sigma

	return Re, Mo, Eo, dx, ut, c_s, dt, nu, mu, rho, delRho, tau
end

function transform(x::Vector)
	x[1] = abs(x[1])
	x[2] = abs(x[2])
	x[3] = abs(x[3])
	return x
end

function targetF(x::Vector)
	Re_t = 25.0
	Eo_t = 5.0
	Mo_t = 1e-4

	x = transform(x)

	g = 		x[1]
	diameter = 	x[2]
	sigma = 	x[3]

	# g = 		abs(x[1])
	# diameter = 	abs(x[2])
	# sigma = 	abs(x[3])

	Re, Mo, Eo = getOtherParams(g,sigma,diameter)
	
	return (Re_t - Re)^2 + (Eo_t - Eo)^2 + (Mo_t - Mo)^2
end

println("functions.jl wurde geladen")