#!/usr/bin/env julia

versioninfo()

####################### FUNCTIONS ######################################

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

out = zeros(size(ARGS,1),2)
header = ["time"  "rho_sum"]

println("num \t rho_sum")

for i = 1 : size(nums,1)
    num = nums[i]
    input_file = getFileNameFromNum(num)
    
    data = readdlm(input_file, ';')[2:end,:]
    
    rho_sum = sum(data[:,4])
    
    println(num," \t ", rho_sum)

    out[i,:] = [num rho_sum]
end

A = [header;out]
writedlm("BubbleMass.csv", A, ';')

println("\n\n##################")