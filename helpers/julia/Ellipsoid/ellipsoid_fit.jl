#!/usr/bin/env julia

versioninfo()

using Optim

function heaviside(x)
    x = x > 0.0 ? x : 0
end

function getM(alpha)
    M = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1]
end

function getD(a,c,r)
    b = r^3/(a*c)
    D = [a^(-2) 0 0; 0 b^(-2) 0; 0 0 c^(-2)]
end

function constructA(x,c)
    M = getM(x[1])
    D = getD(x[2],c,15.0)
    A = M * D * inv(M)
    return A
end

function getCenterPos(data)
    x_m = [0.0,0.0,0.0]
    rho_ges = 0
    for i = 1 : size(data,1)
        x_m = x_m + data[i,1:3] * data[i,4]
        rho_ges += data[i,4]
    end
    return x_m /rho_ges
end

function centerData(data)
    v = getCenterPos(data)
    data_new = zeros(size(data))
    for i = 1 : size(data,1)
        data_new[i,1:3] = data[i,1:3] - v
        data_new[i,4] = data[i,4]
    end
    return data_new
end

function objectiveFunction(x,data,c)
    A = constructA(x,c)
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

function StringToInt(s)
    s = reverse(s)
    sum = 0
    for i = 1 : length(s)
        sum += CharToInt(s[i]) * 10^(i-1)
    end
    return sum
end


out = zeros(size(ARGS,1),9)
header = ["time" "alpha_1" "a_1" "b_1" "c_1" "alpha_2" "a_2" "b_2" "c_2"]

for i = 1 : size(ARGS,1)
    input_file = ARGS[i]
    num = input_file[11:(end-4)]
    num = StringToInt(num)
    println("\n",input_file)

    const alpha = 0.0
    const a = 15.0
    const initial = [alpha, a]
    const r = 15.0
    const c = 15.0

    data = centerData(readdlm(input_file, ';')[2:end,:])

    res_1 = optimize(x -> objectiveFunction(x,data,c), initial)
    show(res_1)
    x_1 = Optim.minimizer(res_1)
    println("\n##################")
    println("alpha_1 = ", x_1[1], " -> ", x_1[1]/(2 * pi)*360,"°", "\na_1 = ",x_1[2])
    println("##################\n")    
    res_2 = optimize(x -> objectiveFunction(x,data,x[2]), initial)
    show(res_2)
    x_2 = Optim.minimizer(res_2)
    println("##################")
    println("alpha_2 = ", x_2[1], " -> ", x_2[1]/(2 * pi)*360,"°", "\na_2 = ",x_2[2])
    println("##################")

    #               alpha_1 a_1    b_1          c_1 alpha_2 a_2     b_2           c_2
    out[i,:] = [num x_1[1] x_1[2] r^3/(x_1[2]*c) c x_2[1] x_2[2] r^3/(x_2[2]^2) x_2[2]]
end

out = sortrows(out, by= x -> (x[1]))# [lt=<comparison>,] [rev=false])

A = [header;out]
writedlm("Ellipsoid.csv", A, ';')

println("\n\n##################")