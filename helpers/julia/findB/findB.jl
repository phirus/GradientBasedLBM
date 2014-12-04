chi = 5;
B = zeros(19)
B[1] = (2+2*chi)/(3*chi +12)
B[2:7] = chi / (6*chi + 24)
B[8:19] = 1 / (6*chi + 24)

show(B)