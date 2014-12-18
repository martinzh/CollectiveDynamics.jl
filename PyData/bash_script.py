
import sys
  
nombre = '../Julia/runs_script.sh'

script = open(nombre,'w')

N = 12
s = 0.05

script.write("#!/bin/bash\n")

for i in range(1,N+1):
  if i%8 == 0:
    script.write( "julia simul.jl " + repr(round(i*s,2)) + "\n")
  else:
    script.write( "nohup julia simul.jl " + repr(round(i*s,2)) + " &\n")


script.close()