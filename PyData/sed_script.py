# Script para modificar las "," por "\t" en las trayectorias y velocidades


import sys

# nombre = '../sed_script.sh'

ruta  = sys.argv[1]
delm  = sys.argv[2]

start = float(sys.argv[3])
stop  = float(sys.argv[4])
step  = float(sys.argv[5])

N     = int((stop - start)/step)

# nombre = '../../'+ ruta +'/sed_script.sh'
nombre = ruta +'/sed_script.sh'
# nombre = '../../DATA/sed_script.sh'

script = open(nombre,'w')


script.write("#!/bin/bash\n")

for i in range(N+1):
  # script.write( "sed 's/,/\t/g' " + "data_f" + repr(round(i*s,2)) + "/trays.txt > " + "data_f" + repr(round(i*s,2)) + "/trays_mod.txt\n")
  # script.write( "sed 's/,/\t/g' " + "data_f" + repr(round(i*s,2)) + "/vels.txt > "  + "data_f" + repr(round(i*s,2)) + "/vels_mod.txt\n")

  script.write( "sed 's/"+delm+"/\t/g' " + "data_f" + repr(round(i*step,6)) + "/trays.txt > " + "data_f" + repr(round(i*step,6)) + "/trays_mod.txt\n")
  script.write( "sed 's/"+delm+"/\t/g' " + "data_f" + repr(round(i*step,6)) + "/vels.txt > "  + "data_f" + repr(round(i*step,6)) + "/vels_mod.txt\n")


script.close()
