
phi_1 = [1 1 1  1  1  1  1  1  1]
phi_2 = [0 1 2  1  2  1  2  1  2]
phi_3 = [0 1 4  1  4  1  4  1  4]
phi_4 = [0 1 1  0 -1 -1 -1  0  1]
phi_5 = [0 1 2  0 -2 -1 -2  0  2]
phi_6 = [0 0 1  1  1  0 -1 -1 -1]
phi_7 = [0 0 2  1  2  0 -2 -1 -2]
phi_8 = [0 1 0 -1  0  1  0 -1  0]
phi_9 = [0 0 1  0 -1  0  1  0 -1]

w = [phi_1 , phi_2 , phi_3 , phi_4 , phi_5 , phi_6 , phi_7 , phi_8 , phi_9]
show(w)
println("\n")

v = zeros(size(w));
show(v)
println("\n")
 
v[1,:] = w[1,:]
show(v)
println("\n") 

for n = 2 : 9

    sum = zeros(w[1,:])

    for i = 1 : (n-1)
        sum = sum + ( [v[i,:] * w[n,:]'] / [v[i,:] * v[i,:]'] * v[i,:] )
    end
    v[n,:] = w[n,:] - sum

    show(v)
    println("\n") 

    if (n == 2)
        v[2,:] = round(v[2,:] * 3)
    end

    if (n == 3)
        v[3,:] = round(v[3,:] * 9/2)
    end

    if (n == 5)
        v[5,:] = round(v[5,:] * 3)
    end

    if (n == 7)
        v[7,:] = round(v[7,:] * 3)
    end
end

 
v_inv = inv(v)
show(round(v_inv * 36))
println("\n") 