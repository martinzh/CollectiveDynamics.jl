# Este script acepta entrada de la linea de comandos
# 1 -> folder en el que se encuentran los datos
# 2 -> valor de la probabilidad
# 3 -> incremento en el valor del campo
# 4 -> campo inicial
# 5 -> campo final
# 6 -> grafica : 1 -> vel 0 -> dr2

import sys

nombre = '../../DATA_1r/DATA/graf.gp'

script = open(nombre,'w')

#script.write("set terminal epslatex size 25.0cm,15cm standalone color colortext 10\n")
#script.write("set output 'test.tex'\n")

script.write("set terminal png size 800,600\n")

# script.write("set output 'clust_"+prob+".png'\n")
script.write("set output 'trayectorias.png'\n")

# script.write("set key top right\n")
# script.write("set xlabel 'tiempo'\n")
script.write("set grid\n")
script.write("unset key\n")

# script.write("set pointsize 0.1\n")
# script.write("set yrange [0:1e5]\n")
#script.write("set label 1 '$p_c$' at 50000,1100 center \n")

# script.write("set xrange [8e8:1e9]\n")
# script.write("set ylabel 'tamano/N'\n")

N = 1080

data = '"data_f'+sys.argv[1]+'/trays_mod.txt"'

script.write('plot ')

for i in range(2,N-2,2):

	script.write(data + " u " + repr(i) + ":" + repr(i+1) + " w lp , ")
# script.write('\n')

script.close()
