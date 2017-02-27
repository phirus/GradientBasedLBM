#!/usr/bin/env julia

versioninfo()

using Optim

include("functions.jl")


out = zeros(size(ARGS,1),3)
header = ["time" "alpha" "a"]

for i = 1 : size(ARGS,1)
    input_file = ARGS[i]
    num = input_file[11:(end-4)]
    num = StringToInt(num)
    println("\n",input_file)

    alpha = 0.0
    a = 15.0
    initial = [alpha, a]

    data = readdlm(input_file, ';')[2:end,:]
    res = optimize(x -> targetF(x,centerData(data)), initial)
    show(res)
    x = Optim.minimizer(res)

    println("\n##################")
    println("alpha = ", x[1], " -> ", x[1]/(2 * pi)*360,"°", "\na = ",x[2])
    println("\n##################")

    out[i,:] = [num x[1] x[2]]
end

out = sortrows(out, by= x -> (x[1]))# [lt=<comparison>,] [rev=false])

A = [header;out]
writedlm("out.csv", A, ';')

println("\n\n##################")
