function heaviside(x)
	x = x > 0.0 ? x : 0
end

function getM(alpha)
	M = [cos(alpha)  sin(alpha) 0;
		 -sin(alpha) cos(alpha) 0;
		           0          0 1]
end

function getD1(a,r)
	b = getB1(a,r)
	D = [a^(-2)      0     0;
		      0 b^(-2)     0;
		      0      0 r^(-2)]
end

function getD2(a,r)
	b = getB2(a,r)
	D = [a^(-2)      0     0;
		      0 b^(-2)     0;
		      0      0 a^(-2)]
end

function getB1(a,r)
	b = r^2/a
end

function getB2(a,r)
	b = r^3/a^2
end

function constructA(alpha,a)
	M = getM(alpha)
	D = getD2(a,15.0)
	A = M * D * inv(M)
	return A
end

function getCenterPos(data)
	x_m = 0
	y_m = 0
	z_m = 0
	rho_ges = 0
	for i = 1 : size(data,1)
		rho = data[i,4]
		x_m += data[i,1] * rho
		y_m += data[i,2] * rho
		z_m += data[i,3] * rho
		rho_ges += rho
	end
	x_m /= rho_ges
	y_m /= rho_ges
	z_m /= rho_ges

	return result = [x_m y_m z_m]
end

function centerData(data)
	v = getCenterPos(data)

	data_new = zeros(size(data))

	for i = 1 : size(data,1)
		data_new[i,1] = data[i,1] - v[1]
		data_new[i,2] = data[i,2] - v[2]
		data_new[i,3] = data[i,3] - v[3]
		data_new[i,4] = data[i,4]
	end

	return data_new
end

function targetF(x,data)
	alpha = x[1]
	a = x[2]
	A = constructA(alpha,a)
	sum = 0
	for i = 1 : size(data,1)
		sum += heaviside((transpose(data[i,1:3]) * A * data[i,1:3])[1])^2 * data[i,4]
	end
	return sum
end

function CharToInt(c)
	if c == '1'
		return 1
	elseif c == '2'
		return 2
	elseif  c == '3'
		return  3
	elseif c == '4'
		return 4
	elseif c == '5'
		return 5
	elseif c == '6'
		return 6
	elseif c == '7'
		return 7
	elseif c == '8'
		return 8
	elseif c == '9'
		return 9
	else
		return 0
	end
end

function tail(s)
	if length(s) <= 1
     	return ""
    else
		return s[2:end]
	end
end

function StringToInt(s)
	if length(s) == 0
		return 0
	else
		return CharToInt(s[1]) * 10^(length(s)-1) + StringToInt(tail(s))
	end
end


println("functions.jl wurde geladen")
