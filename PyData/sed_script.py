
import sys

nombre = '../../DATA/sed_script.sh'

script = open(nombre,'w')

N = 20
s = 0.05

script.write("#!/bin/bash\n")

for i in range(1,N+1):
  if i%6 == 0:
    script.write( "julia simul.jl " + repr(round(i*s,2)) + "\n")
  else:
    script.write( "nohup julia simul.jl " + repr(round(i*s,2)) + " &\n")


script.close()
