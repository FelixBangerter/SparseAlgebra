import Base.show
export BCSRSpM

type BCSRSpM
	val::Array{Array{Float64,2},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	pattern::Array{Int,1}
end

function show(io::IO, SpM::BCSRSpM)
	M = fillSpMWithChar(SpM,*)
	display(M)
end

function fillSpMWithChar(SpM::BCSRSpM,a::Any)
	dim::Int32		= SpM.pattern[end]
	M::Array{Any,2}	= fill(a,dim,dim)
	len::Int32		= length(SpM.rowptr)-1
	δ::Array{Int,1}	= [SpM.rowptr[i+1]-SpM.rowptr[i] for i in 1:len]
	
	j = 1
	for i = 1:len
		z = 0
		while z < δ[i]
			M[pattern[2*(i-1)+1]:(pattern[2*(i-1)+2]),pattern[2*(SpM.col[j]-1)+1]:(pattern[2*(SpM.col[j]-1)+2])] = SpM.val[j]
			j += 1
			z += 1
		end
	end

	return M
end

function computeBlock(M::Array{Float64,2},pattern::Array{Int,1},i::Int,j::Int)
	block = M[pattern[2*(i-1)+1]:(pattern[2*(i-1)+2]),pattern[2*(j-1)+1]:(pattern[2*(j-1)+2])]
	return block
end 

function convertMToSpMBCSR(M::Array{Float64,2},pattern::Array{Int,1})
	
	SpM::BCSRSpM	= BCSRSpM([],[],[0],pattern)
	s::Int32		= length(pattern)/2 
	nnzb			= 0
	
	for i = 1:s
		for j = 1:s
			block = computeBlock(M,pattern,i,j)
			if norm(block) != 0.
				push!(SpM.val, block)
				push!(SpM.col,j)
				nnzb += 1
			end
		end
		push!(SpM.rowptr,nnzb)
	end

	return SpM
end

function convertSpMToMBCSR(SpM::BCSRSpM)
	M = fillSpMWithChar(SpM,0.)
	return M
end
