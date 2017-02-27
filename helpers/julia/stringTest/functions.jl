
function CharToInt(c)
	if c == '1'
		return 1
	elseif c == '2'
		return 2
	elseif  c == '3'
		return  3
	elseif c == '4'
		return 4
	elseif c == '5'
		return 5
	elseif c == '6'
		return 6
	elseif c == '7'
		return 7
	elseif c == '8'
		return 8
	elseif c == '9'
		return 9
	else
		return 0
	end
end

function tail(s)
	if length(s) <= 1
		return ""
	else
		return s[2:end]
	end
end

function StringToInt(s)
	if length(s) == 0
		return 0
	else
		return CharToInt(s[1]) * 10^(length(s)-1) + StringToInt(tail(s))
	end
end

println("functions.jl wurde geladen")
