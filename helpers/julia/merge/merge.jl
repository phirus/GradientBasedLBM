#!/usr/bin/env julia

input_file_1 = ARGS[1]
data_1 = readdlm(input_file_1, ';')
input_file_2 = ARGS[2]
data_2 = readdlm(input_file_2, ';')

A = [data_1 data_2[:,2:end] ]
writedlm("out.csv", A, ';')