import Base.show
import Base.*
import QuantumLab
export BCSRSpM

#Sparse matrix type in block compressed sparse row format
type BCSRSpM
	val::Array{Array{Float64,2},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	rowpattern::Array{Int,1}
	colpattern::Array{Int,1}
end

#Sparse matrix type in block compressed sparse row format for symmetric matrices
type symBCSRSpM
	val::Array{Union{Array{Float64,2},LowerTriangular{Float64,Array{Float64,2}}},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	pattern::Array{Int,1}
end

*(SpM::BCSRSpM,vec::Array{Float64,1}) = multiplySpMV(SpM,vec)

#Funktion generiert Pattern aus Basis nach Shells
function computePatternFromBasis(basis::Array{QuantumLab.ShellModule.Shell,1})
    nbf     = []
    pattern = []
    sum     = 1
    [push!(nbf,QuantumLab.ShellModule.nbf(basis[i])) for i in 1:length(basis)]
    δ = [nbf[i+1]-nbf[i] for i in 1:length(nbf)-1]
    δ = vcat(0,δ)
    [if(δ[i]<0) δ[i] = 0. end for i in 1:length(nbf)]
    
    for i = 1:length(nbf)
        push!(pattern,sum)
        push!(pattern,sum+δ[i])
        sum += nbf[i]
    end

    return pattern
end

function show(io::IO, SpM::BCSRSpM)
	M = fillSpMWithChar(SpM,*)
	display(M)
end

#Converts a LowerTriangular to the corresponding symmetric matrix
function symmetrizeLowerTriangular(lt::LowerTriangular{Float64,Array{Float64,2}})
	M::Array{Float64,2} = Array{Float64,2}(lt)
	n					= size(M,1)
	N					= zeros(n,n)
	for i = 1:n, j = 1:n
		N[i,j] = M[i,j]
		if i == j N[i,j] = 0. end
	end

	return M+N'
end

#Computes the dimension of block to be stored
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

function convertSpMToMBCSR(SpM::BCSRSpM)
	M = fillSpMWithChar(SpM,0.)
	return M
end

function computeSegmentForBlock(SpM::Union{BCSRSpM,symBCSRSpM},i,j)
	if isa(SpM,BCSRSpM)
		rowpattern = SpM.rowpattern
		colpattern = SpM.colpattern
	elseif isa(SpM,symBCSRSpM)
		rowpattern = SpM.pattern
		colpattern = SpM.pattern
	end
	a = rowpattern[2*(i-1)+1]
	b = rowpattern[2*(i-1)+2]
	c = colpattern[2*(SpM.col[j]-1)+1]
	d = colpattern[2*(SpM.col[j]-1)+2]

	return a,b,c,d
end
	
function convertSpMToSMBCSR(SpM::symBCSRSpM)
	dim::Int32			= SpM.pattern[end]
	M::Array{Any,2}		= fill(*,dim,dim)
	len::Int32			= length(SpM.rowptr)-1
	δ::Array{Int,1}		= [SpM.rowptr[i+1]-SpM.rowptr[i] for i in 1:len]

	j = 1
	for i = 1:len
		z = 0
		while z < δ[i]
			a,b,c,d = computeSegmentForBlock(SpM,i,j)
			el = SpM.val[j]
			if a == c && b == d && typeof(SpM.val[j]) == LowerTriangular{Float64,Array{Float64,2}} el = symmetrizeLowerTriangular(el) end
			M[a:b,c:d] = el
			M[c:d,a:b] = el
			j += 1
			z += 1
		end
	end

	return M			
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
			a,b,c,d = computeSegmentForBlock(SpM,i,j)
			M[a:b,c:d] = SpM.val[j]
			j += 1
			z += 1
		end
	end

	return M
end

function computeDifferenceRowptr(i::Int,rowptr)
    a = rowptr[i]+1
    b = rowptr[i+1]
	return a,b
end

#Apply pattern of sparse matrix to a vector
function convertVToBV(vec::Array{Float64,1},SpM::BCSRSpM)
	len::Int64							= Int64(round(length(SpM.colpattern)/2,0))
	blockvec::Array{Array{Float64,1}}	= []
	for i = 1:len
		push!(blockvec,vec[SpM.colpattern[2*(i-1)+1]:SpM.colpattern[2*(i-1)+2]])
	end
		
	return blockvec
end

function multiplySpMV(SpM::BCSRSpM,vec::Array{Float64,1})
	len						= length(SpM.rowptr)-1
	res::Array{Float64,1}	= zeros(length(vec)) 
	blockVec	= convertVToBV(vec,SpM)
	
	for i = 1:len
		a,b = computeDifferenceRowptr(i,SpM.rowptr)
		for j = a:b
			res[SpM.colpattern[2*(i-1)+1]:SpM.colpattern[2*(i-1)+2]] += SpM.val[j] * blockVec[SpM.col[j]]
		end
	end
	
	return res
end

#function multiplySpMM(SpM::BCSRSpM,M::Array{Float64,2})
#	len = length(SpM.rowptr) -1
#	current_row = 1
#	res = zeros(M)
#
#	for i = 1:len
#		a,b = computeDifferenceRowptr(i,SpM.rowptr)
#		for j = a:b
#			for k = 1:len
#				res[k,col[j]] += val[j] * M[k,current_row]
#			end
#		end
#		current_row += 1
#		
#		return res
#end
