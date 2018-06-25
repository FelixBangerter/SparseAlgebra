function CSR(M::Array{Float64,2})
	len::Int32				= length(M)
	s::Int32				= size(M,1)
	val::Array{Float64,1}	= []
	col::Array{Int32,1}		= [] 
	rowptr::Array{Int32,1}	= [0]

	vecM					= convertMatrixToVector(M,s)
	val,col,rowptr			= convertDenseToSparse(len,vecM,s,val,col,rowptr,M)
	
	println("Val",val)
	println("Col",col)
	println("RowPtr",rowptr)
end

function convertDenseToSparse(len::Int32,vecM,s,val::Array{Float64,1},col::Array{Int32,1},rowptr::Array{Int32,1},M::Array{Float64,2})
	j = 0
	for i = 1:len
		if vecM[i] != 0
			push!(val,vecM[i])
			push!(col,i-j*s)
		end
		if i%s == 0	j += 1 
			push!(rowptr,rowptr[j]+countnz(M[(j):(j), :]))
		end
	end
	return val,col,rowptr
end
	
function convertMatrixToVector(M::Array{Float64,2},s::Int32)
	vecM::Array{Float64,1} = []
	for i = 1:s, j = 1:s
		push!(vecM,M[i,j])
	end
	return vecM
end
