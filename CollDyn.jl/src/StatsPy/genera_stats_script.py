# Script para automatizar estadisticas en todas las carpetas


import sys

nombre = './stats_script.sh'
# nombre = '../sed_script.sh'

script = open(nombre,'w')

start = float(sys.argv[1])
stop  = float(sys.argv[2])
step  = float(sys.argv[3])
N     = int((stop - start)/step)

script.write("#!/bin/bash\n")

for i in range(1,N+1):
  script.write( "python Stats.py " + repr(round(i*step,4)) + "\n")

script.close()
