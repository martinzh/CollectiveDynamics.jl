## =========================== ## ## =========================== ##

## 				Libreria para estadistica de datos de simulacion			 ##	
##										 Martin Zumaya Hernandez 									 ##
## 														9/2/15														 ##

## =========================== ## ## =========================== ##

module Tmp
using DataFrames
using Wavelets

## =========================== ## ## =========================== ##
# 										  		INICIA MODULO
## =========================== ## ## =========================== ##


# ext = ".dat" #Extension de los archivos de datos
ext = ".txt" #Extension de los archivos de datos
dists = "dists/" #Directorio con archivos de distancias

## =========================== ## ## =========================== ##

# Obtiene parametros de archivo en el directorio path

function GetParams(path::String)

	data = "$path/params$ext"

	println(data)

	raw_params = readdlm(data,'=')

	# println(typeof(raw_params))
	
	# println(raw_params)

	# for i = 1:size(raw_params,1)
	# 	println("$i $(raw_params[i,:][2])")
	# end

	# println(convert(Array{Float64,1}, raw_params[1:end,2]))
	return convert(Array{Float64,1}, raw_params[1:end,2]) 

end

## =========================== ## ## =========================== ##


	# Obtiene trayectorias de archivo

function	GetVecs(path::String, comp::String)

	data = "$path/$comp$ext"

	println(data)

	return readdlm(data)

end

## =========================== ## ## =========================== ##


function	GetDists(path::String, t::Int64)

	# Obtiene matriz de distancias

	dist_folder = "dist_mat"

	data = "$path/$dist_folder/$t$ext"

	println(data)

	raw_dists = readdlm(data,'\t')

	# Matriz traignular de distancias

	dist_LT = tril(convert(Array{Float64,2},raw_dists[1:end,1:end-1]))

	# Vector con de distancias

	vec_dists = vcat(dist_LT[2:end,1]) # primero a mano

	for i = 2:size(dist_LT,1)

		vec_dists = vcat(vec_dists,dist_LT[i+1:end,i])

	end

	return vec_dists

end

## =========================== ## ## =========================== ##

function VecProm(vecs::Array{Float64,2}, itTot::Int64, step::Int64)

	# tiempo = vcat([1],[ i*step for i = 1:itTot])

	tiempo = size(vecs,1)
	NumCols = size(vecs,2)

	# println("$tiempo,$NumCols")

	# VecsProm = Array(Array{Float64,1}, tiempo)

	VecsProm = zeros(tiempo,2)

	for i in 1:tiempo
		vProm = zeros(2)
		# println(vProm)
		# k = 0
		for j in 1:2:NumCols
			
			# vProm += [vecs[i,j],vecs[i,j+1]]

			vProm[1] += vecs[i,j]
			vProm[2] += vecs[i,j+1]

			# println("$i\t$j\t$(j+1)")
			# broadcast!(+,vProm,vProm,[vecs[i,j],vecs[i,j+1]])
			# println(vProm)
			# k += 1
		end
		# println(k)
		# println(0.5*NumCols)

		scale!(vProm,1.0/(0.5*NumCols))
		# println(vProm)
		
		# push!(VecsProm, vProm)
		
		VecsProm[i,1] = vProm[1]
		VecsProm[i,2] = vProm[2]

	end

	return VecsProm
end

## =========================== ## ## =========================== ##

function OrderParam(vecs::Array{Float64,2}, v0::Float64)

	tiempo = size(vecs,1)
	NumParts = size(vecs,2)

	psi = zeros(tiempo)

	for i in 1:tiempo

		vProm = zeros(2)

		for j in 1:2:NumParts
			vProm += [vecs[i,j],vecs[i,j+1]]
			# broadcast!(+,vProm,vProm,[vecs[i,j],vecs[i,j+1]])
		end

		psi[i] = norm(vProm)/(0.5*NumParts*v0)

	end

	return psi
end

## =========================== ## ## =========================== ##
#								 					TERMINA MODULO 
## =========================== ## ## =========================== ##

end