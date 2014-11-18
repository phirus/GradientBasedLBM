
phi_Humi_1  = [  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]
phi_Humi_2  = [-30 -11 -11 -11 -11 -11 -11   8   8   8   8   8   8   8   8   8   8   8   8]
phi_Humi_3  = [ 12  -4  -4  -4  -4  -4  -4   1   1   1   1   1   1   1   1   1   1   1   1]
phi_Humi_4  = [  0   1  -1   0   0   0   0   1  -1   1  -1   1  -1   1  -1   0   0   0   0]
phi_Humi_5  = [  0  -4   4   0   0   0   0   1  -1   1  -1   1  -1   1  -1   0   0   0   0]
phi_Humi_6  = [  0   0   0   1  -1   0   0   1   1  -1  -1   0   0   0   0   1  -1   1  -1]
phi_Humi_7  = [  0   0   0  -4   4   0   0   1   1  -1  -1   0   0   0   0   1  -1   1  -1]
phi_Humi_8  = [  0   0   0   0   0   1  -1   0   0   0   0   1   1  -1  -1   1   1  -1  -1]
phi_Humi_9  = [  0   0   0   0   0  -4   4   0   0   0   0   1   1  -1  -1   1   1  -1  -1]
phi_Humi_10 = [  0   2   2  -1  -1  -1  -1   1   1   1   1   1   1   1   1  -2  -2  -2  -2] 
phi_Humi_11 = [  0  -4  -4   2   2   2   2   1   1   1   1   1   1   1   1  -2  -2  -2  -2]
phi_Humi_12 = [  0   0   0   1   1  -1  -1   1   1   1   1  -1  -1  -1  -1   0   0   0   0]
phi_Humi_13 = [  0   0   0  -2  -2   2   2   1   1   1   1  -1  -1  -1  -1   0   0   0   0]
phi_Humi_14 = [  0   0   0   0   0   0   0   1  -1  -1   1   0   0   0   0   0   0   0   0]
phi_Humi_15 = [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1  -1  -1   1]
phi_Humi_16 = [  0   0   0   0   0   0   0   0   0   0   0   1  -1  -1   1   0   0   0   0]
phi_Humi_17 = [  0   0   0   0   0   0   0   1  -1   1  -1  -1   1  -1   1   0   0   0   0]
phi_Humi_18 = [  0   0   0   0   0   0   0  -1  -1   1   1   0   0   0   0   1  -1   1  -1]
phi_Humi_19 = [  0   0   0   0   0   0   0   0   0   0   0   1   1  -1  -1  -1  -1   1   1]

Humi1 = [phi_Humi_1 ,phi_Humi_2 ,phi_Humi_3 ,phi_Humi_4 ,phi_Humi_5 ,phi_Humi_6 ,phi_Humi_7 ,phi_Humi_8 ,phi_Humi_9 ,phi_Humi_10,phi_Humi_11,phi_Humi_12,phi_Humi_13,phi_Humi_14,phi_Humi_15,phi_Humi_16,phi_Humi_17,phi_Humi_18,phi_Humi_19]

show(Humi1)
println("\n")

# adpat to ohter ordering of velocities
vergleich = [Humi1[:,1] Humi1[:,2] Humi1[:,8] Humi1[:,4] Humi1[:,9] Humi1[:,3] Humi1[:,11] Humi1[:,5] Humi1[:,10] Humi1[:,6] Humi1[:,12] Humi1[:,16] Humi1[:,13] Humi1[:,17] Humi1[:,7] Humi1[:,14] Humi1[:,18] Humi1[:,15] Humi1[:,19] ]
show(vergleich)
println("\n")
