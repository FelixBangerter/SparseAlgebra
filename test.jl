import QuantumLab

#Funktion generiert Pattern aus Basis nach Shells
function test(basis::Array{QuantumLab.ShellModule.Shell,1})
	nbf		= []
	pattern	= []
	sum		= 1
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
