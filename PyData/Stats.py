# Manejo de los datos de la simulacion con Python
# Libreria para automatizar
# Martin Zumaya Hernandez 2014

from LibData import *
import pylab
import sys
# import numpy as np

# ruta a directorio principal, se identifica por f

# f = 0.2 #conectividad en funcion de fraccion de N
f = float(sys.argv[1]) #conectividad en funcion de fraccion de N

#path = "../data_f"+ str(f) +"/"
#path = "../DATA/data_f"+ str(f)+ "_ro"+str(p) +"/"
# path = "../GitRepos/DATA/data_f"+ str(f) +"/"

path = "../../DATA/data_f"+ str(f) +"/"

# img_path = "../../Grafs/con_r0/datos_clusters"
img_path = "../../Grafs/con_rprom/"
# img_path = "../../Grafs/"+ str(f) +"/"

# print (path)

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
	# paso

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
step     = params[9] #frecuencia de muestreo

# step    = 250 #frecuencia de muestreo (Se obtiene del archivo de params)
num_bin = 200

bin_vec = arange(1,num_bin+1)

tiempo = arange(1,t_f+step,step)

# print(tiempo)

it_tot = t_f/step #muestras totales
# print(it_tot)

max_cls = [] # Para guardar clsuter mas grande

# cls_data = open(img_path+"cls_data"+str(f)+".dat",'w')
cls_data = open(img_path+"/datos_clusters/cls_data"+str(f)+".txt",'w')

for i in range(-1,int(it_tot)):

	if i == -1:
		j = 1
	else:
		j =	(i+1)*int(step)

	# print(repr(j))

	dists = GetDists(path,j) #obtiene distancias
	# print(dists)

	hist = CalcHist(dists,num_bin) #calcula histograma

	print(str(f) + "\t" + repr(j) + "\tr_prom:" + repr(hist[1]) + "\tr_0:" + repr(r_0) + "\tr_max:" + repr(hist[2]))

	# adj = CalcAdjs(dists,r_0) #calcula adjacencias con r_0
	# adj = CalcAdjs(dists,hist[1]) #calcula adjacencias con r_prom
	adj = CalcAdjs(dists,hist[2]) #calcula adjacencias con r_max
	# PrintAdjs(adj)

	# print(sum(hist[0]))
	# print(hist[1:4])

	clusters = Cluster(adj)
	# print(clusters[1])

	# print("mas grande:")
	# print(clusters[2]/N)

	max_cls.append(clusters[2]/N)

	pylab.plot(hist[4]*bin_vec,hist[0])

	cls_data.write(repr(j) + "\t" + repr(clusters[2]/N) + "\n")

pylab.savefig(img_path+"/graficas_distribucion/dist_"+str(f)+".png")

cls_data.close()

# print(max_cls)
# pylab.plot(tiempo,max_cls)
# pylab.savefig(img_path+"cluster_t_"+str(f)+".png")
