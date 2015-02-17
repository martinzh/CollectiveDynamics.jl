
import sys

# nombre = '../Julia/runs_script.sh'
nombre = '../Julia_Obj/runs_script.sh'

script = open(nombre,'w')

start = float(sys.argv[1])
stop  = float(sys.argv[2])
step  = float(sys.argv[3])
T     = sys.argv[4]
rate  = sys.argv[5]
procs = float(sys.argv[6])

N     = int((stop - start)/step)

print(N)

script.write("#!/bin/bash\n")

# for i in range(1,N+1):
# for i in range(N):
#   # if i%6 == 0:
#   if i%(N/2) == 0:
#     script.write( "julia simul.jl " + repr(start + round(i*step,4)) + "\n")
#   else:
#     script.write( "nohup julia simul.jl " + repr(start + round(i*step,4)) + " &\n")

# for i in range(N):
#   # if i%6 == 0:
#   if i%(procs) == 0:
#     script.write( "julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + "\n")
#   else:
#     script.write( "nohup julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + " &\n")

for i in range(N):
  # if i%6 == 0:
  if i%(procs) == 0:
    script.write( "(time julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + ") 2>> tiempo.dat \n")
  else:
    script.write( "(time nohup julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + " &) 2>> tiempo.dat\n")


script.close()
