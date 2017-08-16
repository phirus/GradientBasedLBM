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

function getRadius(data)
    data_new = zeros(size(data,1),5)
    
    for i = 1 : size(data,1)
        data_new[i,1:4] = data[i,1:4] # copy
        data_new[i,5] = sqrt(data[i,1]^2 + data[i,2]^2 + data[i,3]^2) 
    end
    return data_new
end

function myFilter(data,r)
    data_new = zeros(size(data))
    j = 0
    delta = 0.5

    for i = 1 : size(data,1)
        r_tmp = data[i,5] 
        if(r_tmp > r-delta && r_tmp <= r+delta )
            j +=1
            data_new[j,:] = data[i,:]
        end
    end
    return data_new[1:j,:]
end

function averadeRho(data)
    rho = 0
    for i = 1 : size(data,1)
        rho += data[i,4]
    end
    return rho/size(data,1)
end

function objectiveFunction(x,data,rho_fix)
    data = myFilter(data,x[1])
    return (averadeRho(data) - rho_fix)^2
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

out = zeros(size(ARGS,1),10)
header = ["time"  "rho_max" "rho_min" "rho_90" "r_90" "rho_10" "r_10" "rho_50" "r_50" "r_50_"]

println("num \t rho_max \t rho_min \t rho_90 \t r_90 \t rho_10 \t r_10 \t rho_50 \t r_50 \t r_50_")

for i = 1 : size(nums,1)
    num = nums[i]
    input_file = getFileNameFromNum(num)
    
    data = centerData(readdlm(input_file, ';')[2:end,:])
    data = getRadius(data)
    rho_min = minimum(data[:,4])
    rho_max = maximum(data[:,4])

    rho_50 = (rho_min + rho_max) / 2
    rho_90 = (rho_max - rho_min) * 0.9 + rho_min
    rho_10 = (rho_max - rho_min) * 0.1 + rho_min

    result = optimize(r -> objectiveFunction(r,data,rho_90), 8, 14,Brent())
    x = Optim.minimizer(result)
    r_90 = x

    result = optimize(x -> objectiveFunction(x,data,rho_10), 8,14,Brent())
    x = Optim.minimizer(result)
    r_10 = x

    result = optimize(x -> objectiveFunction(x,data,rho_50), 8,14,Brent())
    x = Optim.minimizer(result)
    r_50 = x

    r_50_ = (r_90 + r_10) / 2

    println(num," \t ", rho_max ," \t ", rho_min ," \t ", rho_90 ," \t ", r_90 ," \t ", rho_10 ," \t ", r_10 ," \t ", rho_50 ," \t ", r_50 ," \t ", r_50_)

    out[i,:] = [num rho_max rho_min rho_90 r_90 rho_10 r_10 rho_50 r_50 r_50_]
end

A = [header;out]
writedlm("Spherical.csv", A, ';')

println("\n\n##################")
