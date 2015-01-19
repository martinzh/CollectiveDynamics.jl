# Script para modificar las "," por "\t" en las trayectorias y velocidades


import sys

nombre = '../../DATA/sed_script.sh'
# nombre = '../sed_script.sh'

script = open(nombre,'w')

N = 20
s = 0.05

script.write("#!/bin/bash\n")

for i in range(1,N+1):
  script.write( "sed 's/,/\t/g' " + "data_f" + repr(round(i*s,2)) + "/trays.txt > " + "data_f" + repr(round(i*s,2)) + "/trays_mod.txt\n")
  script.write( "sed 's/,/\t/g' " + "data_f" + repr(round(i*s,2)) + "/vels.txt > "  + "data_f" + repr(round(i*s,2)) + "/vels_mod.txt\n")


script.close()
