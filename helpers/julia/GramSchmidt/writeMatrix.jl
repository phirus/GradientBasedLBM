function createMatrixFile(matrix)

	stream = open("MatrixFile",true,true,true,true,false)

	ymax,xmax = size(matrix)

	for i = 1 : xmax

		for j = 1 : ymax

			write(stream,"matrix[")
			write(stream, string(i-1))
			write(stream,"][")
			write(stream, string(j-1))
			write(stream,"]  = ")
			write(stream, string(matrix[i,j]))
			write(stream, " * t")
			write(stream, " ; \n")
		end
		write(stream, "\n")
	end
	close(stream)
end

println("createMatrixFile wurde geladen")