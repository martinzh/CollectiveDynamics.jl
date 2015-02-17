# Este script acepta entrada de la linea de comandos
# 1 -> folder en el que se encuentran los datos
# 2 -> valor de la probabilidad
# 3 -> incremento en el valor del campo
# 4 -> campo inicial
# 5 -> campo final
# 6 -> grafica : 1 -> vel 0 -> dr2

import sys

nombre = 'graf.gp'

script = open(nombre,'w')

prob = sys.argv[1]

#script.write("set terminal epslatex size 25.0cm,15cm standalone color colortext 10\n")
#script.write("set output 'test.tex'\n")

script.write("set terminal png size 800,600\n")
# script.write("set output 'test.png'\n")

script.write("set output 'clust_"+prob+".png'\n")

script.write("set key top right\n")
script.write("set xlabel 'tiempo'\n")
script.write("set grid\n")
# script.write("set pointsize 0.1\n")
# script.write("set yrange [0:1e5]\n")
#script.write("set label 1 '$p_c$' at 50000,1100 center \n")


# script.write("set xrange [8e8:1e9]\n")
script.write("set ylabel 'tamano/N'\n")

# dp = float(sys.argv[3])
# pi = float(sys.argv[4])
# pf = float(sys.argv[5])
#
# N = int((pf-pi)/dp)

# script.write("fit ")

# for i in range(N+1):
#       e = str(round(abs(pf),2)).rstrip('0').rstrip('.')
#       data = "../"+folder+"/"+prob+"/"+e+"vel.dat"
#       if i != N:
#           script.write('"'+ "a"+repr(i) + "x" + data + '" u 1:(1/$2) w l t '+'"'+e+'"'+", \\\n")
#       else:
#           script.write('"'+ data + '" u 1:(1/$2) w l t '+'"'+e+'"'+"\n")
#       pf -= dp

# script.write('\n')

script.write('plot "cls_data'+ prob +'.dat" w lp')

script.close()
