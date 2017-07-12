#!/usr/bin/env julia

versioninfo()

using Optim

####################### FUNCTIONS ######################################

function heaviside(x)
    x = x > 0.0 ? x : 0
end

function getM(alpha_x,alpha_y,alpha_z)
    
    M_x = [1            0             0;
           0 cos(alpha_x) -sin(alpha_x);
           0 sin(alpha_x)  cos(alpha_x)]

    M_y = [cos(alpha_y) 0 sin(alpha_y);
                      0 1            0;
          -sin(alpha_y) 0 cos(alpha_y)]

    M_z = [cos(alpha_z) -sin(alpha_z) 0; 
           sin(alpha_z)  cos(alpha_z) 0; 
                      0             0 1]
    return M_total = M_x * M_y * M_z
end

function getD(a,b,r)
    c = getC(a,b,r)
    D = [a^(-2) 0 0; 0 b^(-2) 0; 0 0 c^(-2)]
end

function getC(a,b,r)
    c = r^3/(a*b)
end

function constructA(x,r)
    M = getM(x[1],x[2],x[3])
    D = getD(x[4],x[5],r)
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

function objectiveFunction(x,data,r)
    A = constructA(x,r)
    sum = 0
    for i = 1 : size(data,1)
        sum += heaviside((transpose(data[i,1:3]) * A * data[i,1:3])[1])^2 * data[i,4]
    end
    return sum
end

function CharToInt(c)::Int64
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

function StringToInt(s)::Int64
    s = reverse(s)
    sum = 0
    for i = 1 : length(s)
        sum += CharToInt(s[i]) * 10^(i-1)
    end
    return sum
end


function getNumsFromFiles(args)
    nums = zeros(Int64,size(args))
    for i = 1 : size(args,1)
    input_file = args[i]
    num = StringToInt(input_file[11:(end-4)])
    nums[i] = num 
    end

    nums = sort(nums)
end

function getFileNameFromNum(num)
    return string("bubbleFit_",num,".csv")
end

############################ SCRIPT ##################################

nums = getNumsFromFiles(ARGS)

out = zeros(size(ARGS,1),7)
header = ["time" "alpha_x" "alpha_y" "alpha_z" "a" "b" "c"]

println("num \t conv \t calls \t alpha \t beta \t gamma \t a \t b \t c")

const r = 12.0

alpha = 0.0
beta = 0.0
gamma = 0.0
a = 12.0
b = 12.0


for i = 1 : size(nums,1)
    num = nums[i]
    input_file = getFileNameFromNum(num)
    
    initial = [alpha, beta, gamma, a, b]    

    data = centerData(readdlm(input_file, ';')[2:end,:])

    result = optimize(x -> objectiveFunction(x,data,r), initial)
    #show(result)
    x = Optim.minimizer(result)

    alpha = x[1]
    beta = x[2]
    gamma = x[3]
    a = x[4]
    b = x[5]

    println(num," \t ",Optim.converged(result)," \t ",Optim.f_calls(result)," \t ", round(x[1]/(2 * pi)*360,2), "° \t ",round(x[2]/(2 * pi)*360,2),"° \t ",round(x[3]/(2 * pi)*360,2),"° \t ",round(x[4],2)," \t ",round(x[5],2)," \t ",round(r^3/(x[4]*x[5]),2))
    
    out[i,:] = [num alpha beta gamma a b r^3/(a*b)]
end

#out = sortrows(out, by= x -> (x[1]))# [lt=<comparison>,] [rev=false])

A = [header;out]
writedlm("Ellipsoid.csv", A, ';')

println("\n\n##################")