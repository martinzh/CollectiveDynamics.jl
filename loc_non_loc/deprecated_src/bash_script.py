
import sys

nombre = './runs_script.sh'

script = open(nombre,'w')

start = float(sys.argv[1])
stop  = float(sys.argv[2])
step  = float(sys.argv[3])
T     = sys.argv[4]
rate  = sys.argv[5]
eta   = sys.argv[6]
w     = sys.argv[7]


N     = int((stop - start)/step)

if len(sys.argv) == 9 :
	procs = float(sys.argv[8])
else:
	procs = N/2

print(N,procs)

script.write("#!/bin/bash\n")

# for i in range(N):
#   # if i%6 == 0:
#   if i%(procs) == 0:
#     script.write( "julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + "\n")
#   else:
#     script.write( "nohup julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + " &\n")

# for i in range(N):
#   # if i%6 == 0:
#   if i%(procs) == 0:
#     script.write( "(time julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + ") 2>> tiempo.dat \n")
#   else:
#     script.write( "(time nohup julia simul_obj.jl " + repr(start + round(i*step,7)) + " " + T + " " + rate + " &) 2>> tiempo.dat\n")

# eta = "0.35"
# w   = "0.1"

for i in range(N):
  # if i%6 == 0:
  if i%(procs) == 0:
    script.write( "(time julia simul_obj.jl " + repr(i) + " " + T + " " + rate + " " + eta + " " + w + ") 2>> tiempo.dat \n")
  else:
    script.write( "(time nohup julia simul_obj.jl " + repr(i) + " " + T + " " + rate + " " + eta + " " + w + " &) 2>> tiempo.dat\n")


script.close()
