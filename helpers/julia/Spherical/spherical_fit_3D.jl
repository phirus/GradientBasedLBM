#!/usr/bin/env julia

versioninfo()

using Optim

####################### FUNCTIONS ######################################

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

function penalty(rho,rho_fix)
    p = (rho - rho_fix)^2
    return p
end

function objectiveFunction(x,data,rho_fix)
    A = x[1]^(-2) * eye(3)
    sum = 0
    for i = 1 : size(data,1)
        sum += (transpose(data[i,1:3]) * A * data[i,1:3])[1]^2 * penalty(data[i,4],rho_fix)
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
header = ["time"  "rho_max" "rho_min" "rho_90" "r_90" "rho_10" "r_10" "rho_50" "r_50" "r_50_"]

println("num \t conv \t calls \t rho_max \t rho_min \t rho_90 \t r_90 \t rho_10 \t r_10 \t rho_50 \t r_50 \t r_50_")

r_90 = 12.0
r_10 = 12.0

for i = 1 : size(nums,1)
    num = nums[i]
    input_file = getFileNameFromNum(num)

    initial = [alpha, beta, gamma, a, b]

    data = centerData(readdlm(input_file, ';')[2:end,:])
    rho_min = minimum(data[:,4])
    rho_max = maximum(data[:,4])

    rho_50 = (rho_min + rho_max) / 2
    rho_90 = (rho_max - rho_min) * 0.9 + rho_min
    rho_10 = (rho_max - rho_min) * 0.1 + rho_min

    result = optimize(x -> objectiveFunction(x,data,rho_90), initial)
    x = Optim.minimizer(result)
    r_90 = x[1]

    result = optimize(x -> objectiveFunction(x,data,rho_10), initial)
    x = Optim.minimizer(result)
    r_10 = x[1]

    result = optimize(x -> objectiveFunction(x,data,rho_50), initial)
    x = Optim.minimizer(result)
    r_50 = x[1]

    r_50_ = (r_90 + r_10) / 2

    println(num," \t ",Optim.converged(result)," \t ",Optim.f_calls(result)," \t ", rho_max ," \t ", rho_min ," \t ", rho_90 ," \t ", r_90 ," \t ", rho_10 ," \t ", r_10 ," \t ", rho_50 ," \t ", r_50 ," \t ", r_50_)

    out[i,:] = [num rho_max rho_min rho_90 r_90 rho_10 r_10 rho_50 r_50 r_50_]
end

#out = sortrows(out, by= x -> (x[1]))# [lt=<comparison>,] [rev=false])

A = [header;out]
writedlm("Ellipsoid.csv", A, ';')

println("\n\n##################")
