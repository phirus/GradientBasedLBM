require("Humieres.jl")
require("writeMatrix.jl")


phi_1  = [1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1]
phi_2  = [0 1  2  1  2  1  2  1  2  1  2  2  2  2  1  2  2  2  2]
phi_3  = [0 1  4  1  4  1  4  1  4  1  4  4  4  4  1  4  4  4  4]
phi_4  = [0 1  1  0 -1 -1 -1  0  1  0  1  0 -1  0  0  1  0 -1  0]
phi_5  = [0 1  2  0 -2 -1 -2  0  2  0  2  0 -2  0  0  2  0 -2  0]
phi_6  = [0 0  1  1  1  0 -1 -1 -1  0  0  1  0 -1  0  0  1  0 -1]
phi_7  = [0 0  2  1  2  0 -2 -1 -2  0  0  2  0 -2  0  0  2  0 -2]
phi_8  = [0 0  0  0  0  0  0  0  0  1  1  1  1  1 -1 -1 -1 -1 -1]
phi_9  = [0 0  0  0  0  0  0  0  0  1  2  2  2  2 -1 -2 -2 -2 -2]
phi_10 = [0 0 -1 -1 -1  0 -1 -1 -1 -1 -1 -2 -1 -2 -1 -1 -2 -1 -2]
phi_11 = [0 0 -2 -1 -2  0 -2 -1 -2 -1 -2 -4 -2 -4 -1 -2 -4 -2 -4]
phi_12 = [0 0  1  1  1  0  1  1  1 -1 -1  0 -1  0 -1 -1  0 -1  0]
phi_13 = [0 0  2  1  2  0  2  1  2 -1 -2  0 -2  0 -1 -2  0 -2  0]
phi_14 = [0 0  1  0 -1  0  1  0 -1  0  0  0  0  0  0  0  0  0  0]
phi_15 = [0 0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  1]
phi_16 = [0 0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  1  0]
phi_17 = [0 0  1  0 -1  0 -1  0  1  0 -1  0  1  0  0 -1  0  1  0]
phi_18 = [0 0 -1  0 -1  0  1  0  1  0  0  1  0 -1  0  0  1  0 -1]
phi_19 = [0 0  0  0  0  0  0  0  0  0  1 -1  1 -1  0 -1  1 -1  1]

w = [phi_1 ,phi_2 ,phi_3 ,phi_4 ,phi_5 ,phi_6 ,phi_7 ,phi_8 ,phi_9 ,phi_10,phi_11,phi_12,phi_13,phi_14,phi_15,phi_16,phi_17,phi_18,phi_19]
show(w)
println("\n")

v = zeros(size(w));
 
v[1,:] = w[1,:]

for n = 2 : 19

    sum = zeros(w[1,:])

    for i = 1 : (n-1)
        sum = sum + ( [v[i,:] * w[n,:]'] / [v[i,:] * v[i,:]'] * v[i,:] )
    end
    v[n,:] = w[n,:] - sum

    if (n == 2)
        v[2,:] = round(v[2,:] * 19)
    end

    if (n == 3)
        v[3,:] = round(v[3,:] * 21/2)
    end

    if (n == 5)
        v[5,:] = round(v[5,:] * 5)
    end

    if (n == 7)
        v[7,:] = round(v[7,:] * 5)
    end

    if (n == 9)
        v[9,:] = round(v[9,:] * 5)
    end

    if (n == 10)
        v[10,:] = round(v[10,:] * 3)
    end

    if (n == 11)
        v[11,:] = round(v[11,:] * 9)
    end

    if (n == 13)
        v[13,:] = round(v[13,:] * 3)
    end
end

show(v)
println("\n") 

# createMatrixFile(v)
 
test = v - vergleich

show(test)
println("\n") 

v_inv = inv(v)
show(round(v_inv *47880))
# createMatrixFile(round(v_inv *47880))
println("\n") 


test_2 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18, 19]

show(test_2)
show(v * test_2)