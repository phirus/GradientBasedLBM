function getTau(Re)
	tau = (RESOLUTION * MACH_MAX * sqrt(3) /  Re) + 0.5
end

function getParams_alter() 

	gamma = 2;
	delRho = 	RHO_L * (1- 1/ gamma)
	tau = 	getTau(REYNOLDS_MAX_INI)
	dx = 1; 
	c_s = 1/sqrt(3); 
	dt = 1;
	nu = c_s^2 * dt *(tau - 0.5)
	mu = RHO_L * nu

	sigma = sqrt((EOTVOS * (tau - 0.5)^4) / (81 * RESOLUTION^2 * MORTON)) * RHO_L;
	g = (EOTVOS * sigma) / ( RHO_L * (1 - 1/gamma) * RESOLUTION^2 ) ;
	
	Mo = g * mu^4 * delRho/(RHO_L^2 * sigma^3)
	Eo = (g*delRho*RESOLUTION^2)/sigma

	return Mo, Eo, dx, c_s, dt, nu, mu, delRho, tau, sigma, g
end


println("functions.jl wurde geladen")