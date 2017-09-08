#!/usr/bin/env julia

versioninfo()

using Optim

####################### FUNCTIONS ######################################

function getCenterPos(data)
    x_m = [0.0,0.0,0.0]
    rho_ges = 0
    for i = 1 : size(data,1)
        x_m = x_m + data[i,1:3] * data[i,5]
        rho_ges += data[i,5]
    end
    return x_m /rho_ges
end

function centerData(data)
    v = getCenterPos(data)
    data_new = zeros(size(data))
    for i = 1 : size(data,1)
        data_new[i,1:3] = data[i,1:3] - v
        data_new[i,4:6] = data[i,4:6]
    end
    return data_new
end

function getRadius(data)
    data_new = zeros(size(data,1),7)
    
    for i = 1 : size(data,1)
        data_new[i,1:6] = data[i,1:6] # copy
        data_new[i,7] = sqrt(data[i,1]^2 + data[i,2]^2 + data[i,3]^2) 
    end
    return data_new
end

function myFilter(data,r)
    data_new = zeros(size(data))
    j = 0
    delta = 0.5

    for i = 1 : size(data,1)
        r_tmp = data[i,7] 
        if(r_tmp > r-delta && r_tmp <= r+delta )
            j +=1
            data_new[j,:] = data[i,:]
        end
    end
    return data_new[1:j,:]
end

function averadePsi(data)
    psi = 0
    for i = 1 : size(data,1)
        psi += data[i,6]
    end
    return psi/size(data,1)
end

function objectiveFunction(x,data,psi_fix)
    data = myFilter(data,x[1])
    return (averadePsi(data) - psi_fix)^2
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

out = zeros(size(ARGS,1),5)
header = ["time"  "r_p" "r_m" "r_0" "r_0_"]

println("num \t r_p \t r_m \t r_0 \t r_0_")

for i = 1 : size(nums,1)
    num = nums[i]
    input_file = getFileNameFromNum(num)
    
    data = centerData(readdlm(input_file, ';')[2:end,:])
    data = getRadius(data)

    psi_high = 0.9
    psi_low = -0.9
    psi_0 = 0

    result = optimize(r -> objectiveFunction(r,data,psi_high), 8, 14,Brent())
    x = Optim.minimizer(result)
    r_p = x

    result = optimize(x -> objectiveFunction(x,data,psi_low), 8,14,Brent())
    x = Optim.minimizer(result)
    r_m = x

    result = optimize(x -> objectiveFunction(x,data,psi_0), 8,14,Brent())
    x = Optim.minimizer(result)
    r_0 = x

    r_0_ = (r_p + r_m) / 2

    println(num," \t ", r_p ," \t ", r_m ," \t ", r_0 ," \t ", r_0_)

    out[i,:] = [num r_p r_m r_0 r_0_]
end

A = [header;out]
writedlm("Spherical.csv", A, ';')

println("\n\n##################")
