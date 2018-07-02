import Base.show
export CSRSpM

type BCSRSpM
	val::Array{Array{Float64,2},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
end

function show(io::IO, SpM::BCSRSpM)
	M = fillAnyMatrixWithSpMElements(SpM,*)	
	display(M)
end

#function show(io::IO, x::Array{Array{Float64,2},1})
#	len = size(val,1)
#	display(reshape([val[k][i] for i in eachindex(val[1]), k in eachindex(val)],(2,2,len)))
#end

function convertMToSpMBCSR(A::Array{Float64,2})
	blocksize						= 2
	SpM::BCSRSpM					= BCSRSpM([],[],[0])
	s::Int32						= size(A,1)/blocksize
	tmp								= []
	nnzb = 0
	
	for i = 1:s
		for j = 1:s
		tmp = A[1+(i-1)*blocksize:i*blocksize,1+(j-1)*blocksize:j*blocksize]
			if norm(tmp) != 0.
				push!(SpM.col,j)
				push!(SpM.val,tmp)
				nnzb += 1
			end
		end
		push!(SpM.rowptr,nnzb)		
	end
	return SpM
end

function convertSpMToMBCSR(SpM::BCSRSpM)
	M = fillAnyMatrixWithSpMElements(SpM,0.0)	
	println(M)
	return M 
end

function fillAnyMatrixWithSpMElements(SpM::BCSRSpM,a::Any)
	dim::Int32			= size(SpM.val[1],1)*(length(SpM.rowptr)-1)
	blocksize			= size(SpM.val[1],1)
	M::Array{Any,2}		= fill(a,dim,dim)
	len::Int32			= length(SpM.rowptr)-1
	δ					= [SpM.rowptr[i+1]-SpM.rowptr[i] for i = 1:len]

	j = 1
	for i = 1:len
		z = 0
		while z < δ[i]
			M[1+(i-1)*blocksize:i*blocksize,1+(SpM.col[j]-1)*blocksize:SpM.col[j]*blocksize] = SpM.val[j]
			j += 1
			z += 1
		end
	end
	
	return M
end
