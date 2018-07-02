import Base.show
export BCSRSpM

type BCSRSpM
	val::Array{Array{Float64,2},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	rowpattern::Array{Int,1}
	colpattern::Array{Int,1}
end

function show(io::IO, SpM::BCSRSpM)
	M = fillSpMWithChar(SpM,*)
	display(M)
end

function fillSpMWithChar(SpM::BCSRSpM,a::Any)
	dim1::Int32		= SpM.rowpattern[end]
	dim2::Int32		= SpM.colpattern[end]
	M::Array{Any,2}	= fill(a,dim1,dim2)
	len::Int32		= length(SpM.rowptr)-1
	δ::Array{Int,1}	= [SpM.rowptr[i+1]-SpM.rowptr[i] for i in 1:len]
	
	j = 1
	for i = 1:len
		z = 0
		while z < δ[i]
			M[rowpattern[2*(i-1)+1]:(rowpattern[2*(i-1)+2]),colpattern[2*(SpM.col[j]-1)+1]:(colpattern[2*(SpM.col[j]-1)+2])] = SpM.val[j]
			j += 1
			z += 1
		end
	end

	return M
end

function computeBlock(M::Array{Float64,2},rowpattern::Array{Int,1},colpattern::Array{Int,1},i::Int,j::Int)
	block = M[rowpattern[2*(i-1)+1]:(rowpattern[2*(i-1)+2]),colpattern[2*(j-1)+1]:(colpattern[2*(j-1)+2])]
	return block
end 

function convertMToSpMBCSR(M::Array{Float64,2},rowpattern::Array{Int,1},colpattern::Array{Int,1})
	
	SpM::BCSRSpM	= BCSRSpM([],[],[0],rowpattern,colpattern)
	s1::Int32		= length(rowpattern)/2 
	s2::Int32		= length(colpattern)/2
	nnzb			= 0
	
	for i = 1:s1
		for j = 1:s2
			block = computeBlock(M,rowpattern,colpattern,i,j)
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

type symBCSRSpM
	val::Array{Union{Array{Float64,2},LowerTriangular{Float64,Array{Float64,2}}},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	pattern::Array{Int,1}
end

function convertSMtoSPM(M::Array{Float64,2},pattern::Array{Int,1},tri::Bool)
	
	SpM::symBCSRSpM		= symBCSRSpM([],[],[0],pattern)
	s::Int32		= length(pattern)/2
	nnzb			= 0
	
	for i = 1:s
		for j = 1:s
			if j <= i
				block = computeBlock(M,pattern,pattern,i,j)
				if tri == true && i == j
					block = LowerTriangular(block)
				end
				if norm(block) != 0.
					push!(SpM.val,block)
					push!(SpM.col,j)
					nnzb += 1
				end
			end
		end
		push!(SpM.rowptr,nnzb)
	end
		
	return SpM
end
	
