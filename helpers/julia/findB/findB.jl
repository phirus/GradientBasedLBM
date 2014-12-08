require("functions.jl")

c_x = [0   1   1   0   -1  -1  -1  0   1   0   1   0   -1  0   0   1   0   -1  0]
c_y = [0   0   1   1   1   0   -1  -1  -1  0   0   1   0   -1  0   0   1   0   -1]
c_z = [0   0   0   0   0   0   0   0   0   1   1   1   1   1   -1  -1  -1  -1  -1]

# 2D 
println("\n2D case:\n")
B_2D = zeros(9)
B_2D[1] = -4/27
B_2D[2] = 2/27
B_2D[3] = 5/108
B_2D[4] = 2/27
B_2D[5] = 5/108
B_2D[6] = 2/27
B_2D[7] = 5/108
B_2D[8] = 2/27
B_2D[9] = 5/108
println("B:")
show(B_2D)
println("\n\nsum(B):",sum(B_2D))
println("\nsum(B * c_x):",prodSum(B_2D,c_x))
println("\nsum(B * c_y):",prodSum(B_2D,c_y))

println("\n3rd sum:")
show(ultraSum(B_2D,c_x,c_y))

println("\n2D case:\n")
B_2D_alter = zeros(9)
B_2D_alter[1] = 3/4
B_2D_alter[2] = 3/4
B_2D_alter[3] = 3/4
B_2D_alter[4] = 3/4
B_2D_alter[5] = 3/4
B_2D_alter[6] = 3/4
B_2D_alter[7] = 3/4
B_2D_alter[8] = 3/4
B_2D_alter[9] = 3/4
println("B:")
show(B_2D_alter)
println("\n\nsum(B):",sum(B_2D_alter))
println("\nsum(B * c_x):",prodSum(B_2D_alter,c_x))
println("\nsum(B * c_y):",prodSum(B_2D_alter,c_y))

println("\n3rd sum:")
show(ultraSum(B_2D_alter,c_x,c_y))

# 3D
println("\n\n3D case:\n")

chi = 2;
B_3D = zeros(19)
B_3D[1] = -(2+2*chi)/(3*chi +12)
B_3D[2] = chi / (6*chi + 24)
B_3D[3] = 1 / (6*chi + 24)
B_3D[4] = chi / (6*chi + 24)
B_3D[5] = 1 / (6*chi + 24)
B_3D[6] = chi / (6*chi + 24)
B_3D[7] = 1 / (6*chi + 24)
B_3D[8] = chi / (6*chi + 24)
B_3D[9] = 1 / (6*chi + 24)
B_3D[10] = chi / (6*chi + 24)
B_3D[11] = 1 / (6*chi + 24)
B_3D[12] = 1 / (6*chi + 24)
B_3D[13] = 1 / (6*chi + 24)
B_3D[14] = 1 / (6*chi + 24)
B_3D[15] = chi / (6*chi + 24)
B_3D[16] = 1 / (6*chi + 24)
B_3D[17] = 1 / (6*chi + 24)
B_3D[18] = 1 / (6*chi + 24)
B_3D[19] = 1 / (6*chi + 24)

println("B:")
show(B_3D)
println("\n\nsum(B):",sum(B_3D))
println("\nsum(B * c_x):",prodSum(B_3D,c_x))
println("\nsum(B * c_y):",prodSum(B_3D,c_y))
println("\nsum(B * c_z):",prodSum(B_3D,c_z))

println("\n3rd sum:")
show(ultraSum3(B_3D,c_x,c_y,c_z))