########################################################
# Libreria para estadistica de datos de simulacion
# Martin Zumaya Hernandez 
# 9/2/15
########################################################

module Tmp
using DataFrames
using Wavelets

########################################################
# 			INICIA MODULO
########################################################

ext = ".dat" #Extension de los archivos de datos
dists = "dists_mat/" #Directorio con archivos de distancias

########################################################

	# Obtiene parametros de archivo en el directorio path

function GetParams(path)

	data = "$path/parametros$ext"

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

########################################################

	# Obtiene trayectorias de archivo

function	GetTrays(path)

	data = "$path/trays$ext"

	println(data)

	raw_trays = readdlm(data,'\t')

	return convert(Array{Float64,2},raw_trays[1:end,1:end-1])

end

########################################################

function	GetVels(path)

	# Obtiene velocidades de archivo

	data = "$path/vels$ext"

	println(data)

	raw_vels = readdlm(data,'\t')

	return convert(Array{Float64,2},raw_vels[1:end,1:end-1])

end

########################################################

function	GetDists(path, t::Int64)

	# Obtiene matriz de distancias

	dist_folder = "dist_mat"

	data = "$path/$dist_folder/$t$ext"

	println(data)

	raw_dists = readdlm(data,'\t')

	# Matriz traignular de distancias

	# dist_LT = triu(convert(Array{Float64,2},raw_dists[1:end,1:end-1]))
	dist_LT = tril(convert(Array{Float64,2},raw_dists[1:end,1:end-1]))

	# Vector con de distancias

	vec_dists = vcat(dist_LT[2:end,1]) # primero a mano

	for i = 2:size(dist_LT,1)

		vec_dists = vcat(vec_dists,dist_LT[i+1:end,i])

	end

	return vec_dists

end

########################################################

########################################################
# 					TERMINA MODULO 
########################################################

end