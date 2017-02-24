#!/bin/bash

N="$1"
k="$2"
w="$3"
T="$4"

for i in {1..15}
do
  nohup julia ../src/work_sim.jl $N $k $w $T $i &
done
