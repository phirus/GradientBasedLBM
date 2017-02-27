#!/usr/bin/env julia

#versioninfo()

using Optim

include("functions.jl")

r = 15.0

# start values
alpha = 0.0
a = 10.0

#for i = i : size(ARGS,1)



#input_file = ARGS[1]
input_file = "bubbleFit_10.csv"

num = input_file[11:(end-4)]
show(num)

# Open the output file for writing
data = readdlm(input_file, ';')[2:end,:]

res = optimize(x -> targetF(x,centerData(data)), [alpha, a])#, NelderMead())
show(res)

#x = res.minimum

#alpha = x[1]
#a = x[2]

println("\n\n##################")

#createPFile(Mo, Eo, c_s, gamma, sigma, g)
