# Script para automatizar estadisticas en todas las carpetas


import sys

nombre = './stats_script.sh'
# nombre = '../sed_script.sh'

script = open(nombre,'w')

N = 20
s = 0.05

script.write("#!/bin/bash\n")

for i in range(1,N+1):
  script.write( "python Stats.py " + repr(round(i*s,2)) + "\n")

script.close()
