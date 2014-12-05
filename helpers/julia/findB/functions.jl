function prodSum(list,fac)
    sum = 0
    s = size(list,1)
    for i = 1 : s
        sum = sum + list[i] * fac[i]
    end
    return sum
end

function ultraSum(list,fac_x,fac_y)
    sum = [[0 0] , [0 0]]
    s = size(list,1)
    for i = 1 : s
        e = [fac_x[i],fac_y[i]]
        sum = sum + list[i] * (e * e')
    end
    return sum
end

function ultraSum3(list,fac_x,fac_y,fac_z)
    sum = [[0 0 0] , [0 0 0] , [0 0 0]]
    s = size(list,1)
    for i = 1 : s
        e = [fac_x[i],fac_y[i],fac_z[i]]
        sum = sum + list[i] * (e * e')
    end
    return sum
end

println("functions.jl wurde geladen")