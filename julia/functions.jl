function getTargetReynolds()
	Re = 50
end


function getTargetEotvos()
	Eo = 10
end

function getTargetMorton()
	Mo = 1e-4
end

function getTau(Reynolds)
	tau = (4*sqrt(3)/Reynolds) + 0.5
end

function getSigma(Eotvos, g, delRho,diameter)
	sigma = g * delRho * diameter^2 / Eotvos
end

function getOtherParams(g, diameter)
	rho = 		1.0
	delRho = 	0.5
	tau = 	getTau(getTargetReynolds())
	sigma = getSigma(getTargetEotvos(), g, delRho, diameter)

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

	Re = 4*sqrt(3)/(tau - 0.5)	
	Mo = g * mu^4 * delRho/(rho^2 * sigma^3)
	Eo = (g*delRho*diameter^2)/sigma

	return Re, Mo, Eo, dx, ut, c_s, dt, nu, mu, rho, delRho, tau
end

# function transform(x::Vector)
# 	x[1] = abs(x[1])
# 	x[2] = abs(x[2])
# 	x[3] = abs(x[3])
# 	return x
# end

# function targetF(x::Vector)
# 	Eo_t = getTargetEotvos()
# 	Mo_t = getTargetMorton()

# 	x = transform(x)

# 	g = 		x[1]
# 	sigma = 	x[2]
# 	diameter = 	x[3]
	

# 	Re, Mo, Eo = getOtherParams(g, sigma, diameter)
	
# 	return  ((Eo_t - Eo)/Eo_t)^2 + ((Mo_t - Mo)/Mo_t)^2 
# end

println("functions.jl wurde geladen")