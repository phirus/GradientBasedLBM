function getTau(Re)
	tau = (RESOLUTION * MACH_MAX * sqrt(3) /  Re) + 0.5
end

function getTerminalRiseVelo(sigma, g)
	tmp = 2.14 * sigma / (RHO_L * DIAMETER) + 0.505 *g * DIAMETER
	if (tmp < 0 ) 
		tmp = 0
	end
	ut = sqrt( tmp )
	return ut
end

function getOtherParams(g, sigma,gamma)

	delRho = 	RHO_L * (1- 1/ gamma)
	tau = 	getTau(REYNOLDS_MAX_INI)
	dx = DIAMETER/RESOLUTION
	ut = getTerminalRiseVelo(sigma, g)
	c_s = ut / MACH_MAX;
	dt = dx / (sqrt(3)*c_s)
	nu = c_s^2 * dt *(tau - 0.5)
	mu = RHO_L * nu
	
	Mo = g * mu^4 * delRho/(RHO_L^2 * sigma^3)
	Eo = (g*delRho*DIAMETER^2)/sigma

	return Mo, Eo, dx, ut, c_s, dt, nu, mu, delRho, tau
end

function transform(x::Vector) # made a function in case i need another transformation
	y =  abs(x)
	return y
end

function targetF(x::Vector)
	x = transform(x)

	g = 		x[1]
	sigma = 	x[2]
	gamma = 	x[3]

	Mo, Eo = getOtherParams(g, sigma, gamma)
	
	return  ((EOTVOS - Eo)/EOTVOS)^2 + ((MORTON - Mo)/MORTON)^2 
end

println("functions.jl wurde geladen")