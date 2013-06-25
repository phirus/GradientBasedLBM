
function getOtherParams(x::Vector)
	g = 		x[1]
	rho = 		x[2]
	delRho = 	x[3]
	diameter = 	x[4]
	sigma = 	x[5]
	tau = 		x[6]

	dx = diameter/40
	ut = sqrt(2.14*sigma /(rho*diameter) + 0.505 *g * diameter)
	c_s = 10*ut;
	dt = dx / (sqrt(3)*c_s)
	nu = c_s^2 * dt *(tau - 0.5)
	mu = rho * nu;
	Re = rho * diameter * ut / mu
	Mo = g * mu^4 * delRho/(rho^2 * sigma^3)
	Eo = (g*delRho*diameter^2)/sigma

	return Re, Mo, Eo, dx, ut, c_s, dt, nu, mu
end

function targetF(x::Vector)
	Re_t = 25
	Eo_t = 20;
	Mo_t = 1e-4
	Re, Mo, Eo = getOtherParams(x)

	return (Re_t - Re)^2 + (Eo_t - Eo)^2 + (Mo_t - Mo)^2
end