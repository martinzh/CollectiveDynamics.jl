# Manejo de los datos de la simulacion con Python
# Libreria para automatizar
# Martin Zumaya Hernandez 2014

from LibData import *

# ruta a directorio principal, se identifica por f

f = 0.01 #conectividad en funcion de fraccion de N

#path = "../data_f"+ str(f) +"/"
#path = "../DATA/data_f"+ str(f)+ "_ro"+str(p) +"/"
# path = "../GitRepos/DATA/data_f"+ str(f) +"/"

path = "../../DATA/data_f"+ str(f) +"/"

print (path)

# Estructura parametros: 
	# Particulas
	# densidad
	# radio
	# f
	# ruido geometrico
	# ruido topologico
	# peso geometrico
	# regimen de velocidad
	# iteraciones
	# paso FALTA IMPLEMENTAR

# obtiene parametros de archivo 
params = GetParams(path)
# print(params)

N        = int(params[0]) #Numero de particulas
ro       = params[1] #Densidad
r_0      = params[2] #Radio interaccion (vel_reg)
f        = params[3] #Fraccion de N en largo alcance
noise_sh = params[4] #Ruido corto
noise_lg = params[5] #Ruido largo
rel_weig = params[6] #Peso Relativo vecindades
reg_vel  = params[7] #Regimen de velocidad
t_f      = params[8] #iteraciones totales
# step     = params[9] #frecuencia de muestreo

step    = 250 #frecuencia de muestreo
num_bin = 150

it_tot = t_f/step #muestras totales
# print(it_tot)

dists = GetDists(path,1) #obtiene distancias
# print(dists)

adj = CalcAdjs(dists,r_0) #calcula adjacencias
# PrintAdjs(adj)

hist = CalcHist(dists,num_bin)
# print(sum(hist[0]))
# print(hist[1:4])

clusters = Cluster(adj)
print(clusters[1])

print("numero de clusters:")
print(len(clusters[1]))

print("indice mas alto:")
print(max(clusters[1]))

print("llaves:")
print(sorted(clusters[1].keys()))

print("valores:")
print(sorted(clusters[1].values()))

print("mas grande:")
print(sorted(clusters[1].values())[-1]/N)



