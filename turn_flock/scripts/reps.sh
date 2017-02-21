#!/bin/bash

N="$1"
k="$2"
T="$3"

for i in {1..5}
do
  nohup julia turn_flock_model.jl $N $k $T $i &
done
