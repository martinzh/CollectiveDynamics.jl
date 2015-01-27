# Manejo de los datos de la simulacion con Python
# Libreria para automatizar
# Martin Zumaya Hernandez 2014

import networkx as nx
import matplotlib.pyplot as plt
from pylab import *
import sys

sys.setrecursionlimit(3000)

ext = ".txt" #Extension de los archivos de datos
dists = "dists/" #Directorio con archivos de distancias

########################################################

# Obtiene parametros de archivo en el directorio path

def GetParams(path):

	params = []

	data = open(path + "params.txt",'rb')

	for line in data.readlines():
	    # print(repr(line.split()[0]) +"\t"+ repr(line.split()[-1]))
	    params.append(float(line.split()[-1]))

	data.close()

	return params

########################################################

# Obtiene matriz de distancias completa
def GetDists(path, t, N):

	DD = [] #Arreglo para guradar distancias

	ruta = path + dists + repr(t) + ext
	# print(ruta)

	data = open(ruta,'rb')

	i = 0

	for line in data.readlines():

		# dd = []

		vals = line.split()

		# print(repr(len(vals))+"\n")
		# print(repr(i)+"\n")

		# for k in range(len(vals)):
		for k in range(i+1,N):

			# dd.append(float(vals[k]))
			DD.append(float(vals[k]))

		# DD.append(dd)

		i += 1

	# print(len(DD))

	return DD

########################################################

# Calcula distribucion de distancias
def CalcHist(dists, num_bin):

	bins = arange(1,num_bin+1) #Arreglo de bins
	hist = zeros(num_bin) #inicializa el histograma

	DD = []

	j = 0

	n = len(dists)

	# print(n)

	for i in range(n):

	    for k in range(n-j):

	        DD.append(dists[i][k])

	    j += 1

	DD.sort() #Ordena de menor a mayor

	#print DD[0]
	#print len(DD)

	epsilon = DD[-1]/num_bin #Tamanio de cada bin en funcion del maximo

	for d in DD:
	    count = d/epsilon
	    if count > 0: hist[int(count-1)] += 1.0 #contador de cada bin

	# Calcula distancia promedio
	r_prom = 0

	# for i in range(int(num_bin)):
	for i in range(int(num_bin)):

		# r_prom += (epsilon*(i)*hist[i])/num_bin
		# r_prom += epsilon*(i+1)*hist[i]
		r_prom += epsilon*i*hist[i]

	r_prom *= 1/num_bin

	#Falta determinar bien el r mas probable

	r_max = np.argmax(hist)*epsilon # distancia mas probable
	std_dev = np.std(DD) #desviacion estandar

	return [hist,r_prom,r_max,std_dev,epsilon,sum(hist)]

########################################################
########################################################

# Calcula distribucion de distancias (version numpy)

def CalcHist1(dists, num_bin):

	dists.sort() #Ordena de menor a mayor

	#print dists[0]
	#print len(dists)

	epsilon = dists[-1]/num_bin #Tamanio de cada bin en funcion del maximo

	# for d in dists:
	#     count = d/epsilon
	#     if count > 0: hist[int(count-1)] += 1.0 #contador de cada bin

	hist = np.histogram(dists,num_bin, density = True)
	# hist = np.histogram(dists,num_bin)

	# Calcula distancia promedio
	r_prom = 0

	# for i in range(int(num_bin)):
	for i in range(int(num_bin)):

		# r_prom += (epsilon*(i)*hist[i])/num_bin
		# r_prom += epsilon*(i+1)*hist[i]
		# r_prom += epsilon*i*hist[0][i]
		r_prom += np.diff(hist[1])[i]*i*hist[0][i]

	# r_prom *= 1/num_bin

	#Falta determinar bien el r mas probable

	r_max = np.argmax(hist[0])*epsilon # distancia mas probable
	std_dev = np.std(dists) #desviacion estandar

	return [hist[0],r_prom,r_max,std_dev,epsilon,sum(hist[0])]

########################################################

# Calcula lista de adjacencia a partir de dists y r
def CalcAdjs(dists,r):

	ADJ = [] #adjacencia toda la red

	for i in range(len(dists)):

		adj = [] #adjacencia particula i

		adj.append(i)

		for j in range(len(dists[i])):

			if dists[i][j] > 0 and dists[i][j] <= r :

				adj.append(j)

		ADJ.append(adj)

	return ADJ

########################################################

########################################################

# Calcula lista de adjacencia a partir de dists y r
def CalcAdjs1(N,dists,r):

	ADJ = [] #adjacencia toda la red

	k = 0

	for i in range(N):

		adj = [] #adjacencia particula i

		adj.append(i)

		for j in range(i+1,N):

			if dists[k] > 0 and dists[k] <= r :

				adj.append(j)

			k += 1

		ADJ.append(adj)

	return ADJ

########################################################

# Imprime lista de adjacencia
def PrintAdjs(adj):
	for i in range(len(adj)):
		print(adj[i])

########################################################

# Etiqueta los elementos de la vecindad

def set_label(U,i,c,label):

	L = len(U[i]) # Accede a la vecindad de la part i

	if L != 0: # Si la vecindad no es vacia

		for j in range(L): # para los vecinos de...

			k = U[i][j] # indice de vecino

			# print (k)

			if label[k] == 0:

				label[k] = label[i]

				set_label(U,k,c,label)

########################################################

# CLUSTERING RECURSIVO

def Cluster(neigh):

	c = 1 # etiqueta de cluster inicial

	for i in range(len(neigh)):

		neigh[i].pop(0)

		#print (neigh[i])

	label = zeros(len(neigh)) # vector de etiquetas
	cl_size = {} # diccionario: llave:cluster, valor:tamanio
	#print (labels)

	for i in range(len(neigh)):

		if label[i] == 0:

			cl_size[c] = 1

			label[i] = c

			# print (label)

			set_label(neigh,i,c,label)

			# print (label)
		else: cl_size[label[i]] += 1

		c += 1

	# numero de clusters
	NC = len(cl_size)

	# indice mas alto
	IM = max(cl_size)

	# llaves
	# sorted(cl_size.keys())
	# valores
	# sorted(cl_size.values())

	# cluster mas grande
	CM = sorted(cl_size.values())[-1]

	# print (label)
	return [label,cl_size,CM]
