import Base.show
export CSRSpM

type CSRSpM
	val::Array{Float64,1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
end

function show(io::IO, x::CSRSpM)
	M = CSRSparseToDense(x.val,x.col,x.rowptr)
	star::Array{Any,2} = fill(*,size(M,1),size(M,1))
	for i = 1:size(M,1), j = 1:size(M,1)
		if M[i,j] != 0.
			star[i,j] = M[i,j]
		end
	end
	display(star)
end

function test(val,col,rowptr)
	SpM::CSRSpM = CSRSpM(val,col,rowptr)
	show(SpM)
end

function computeDifferenceRowptr(i::Int,rowptr)
	a = rowptr[i]+1
	b = rowptr[i+1]
	return a,b
end

function purgeSparseMatrix(SpM::CSRSpM)
	Θ = 1e-5
	len = length(SpM.rowptr)-1
	del::Array{Int} = []
	delrow::Array{Int} = []	
	
	for i = 1:len
		a,b = computeDifferenceRowptr(i,SpM.rowptr)
		for j = a:b
			if SpM.val[j] < Θ
				push!(del,j)
				push!(delrow,i)
			end
		end
	end

	deleteat!(SpM.val,del)
	deleteat!(SpM.col,del)
	for i = 1:length(delrow)
		SpM.rowptr[(delrow[i]+1):end] -= 1
	end
	return SpM
end

function CSRDenseToSparse(M::Array{Float64,2})
	len::Int32				= length(M)
	s::Int32				= size(M,1)
	val::Array{Float64,1}	= []
	col::Array{Int,1}		= [] 
	rowptr::Array{Int,1}	= [0]

	vecM					= convertMatrixToVector(M,s)
	val,col,rowptr			= convertDenseToSparse(len,vecM,s,val,col,rowptr,M)

	return val,col,rowptr
end

function convertDenseToSparse(len::Int32,vecM,s,val::Array{Float64,1},col::Array{Int,1},rowptr::Array{Int,1},M::Array{Float64,2})
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

function CSRSparseToDense(val::Array{Float64,1},col::Array{Int,1},rowptr::Array{Int,1})
	len		= length(rowptr)-1
	M		= zeros(len,len)
	δ		= [rowptr[i+1]-rowptr[i] for i = 1:len]
	
	j = 1
	for i = 1:len
		z = 0
		while(z < δ[i])
			M[i, col[j]] = val[j]
			z += 1
			j += 1
		end
	end

	return M
end

function multiplySpMV(val::Array{Float64,1},col::Array{Int,1},rowptr::Array{Int,1},vec::Array{Float64,1})
	len		= length(rowptr)-1
	res		= zeros(length(vec))
	
	for i = 1:len
		a,b = computeDifferenceRowptr(i,rowptr)
		for j = a:b
			res[i] += val[j] * vec[col[j]]
		end
	end

	return res
end

function multiplyMSpM(val::Array{Float64,1},col::Array{Int,1},rowptr::Array{Int,1},M::Array{Float64,2})
	len			= length(rowptr)-1
	current_row = 1
	res			= zeros(len,len)	
	
	for i = 1:len
		for j = rowptr[i]+1:rowptr[i+1]
			for k = 1:len
				res[k,col[j]] += val[j] * M[k,current_row]
			end
		end
		current_row += 1
	end
	
	return res
end
