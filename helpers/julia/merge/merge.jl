#!/usr/bin/env julia

input_file_1 = "BubblePlot.csv"#ARGS[1]
data_1 = readdlm(input_file_1, ';')
input_file_2 = "Ellipsoid.csv"#ARGS[2]
data_2 = readdlm(input_file_2, ';')

A = [data_1 data_2[:,2:end] ]
writedlm("BubblePlot_complete.csv", A, ';')