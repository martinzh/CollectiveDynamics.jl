#!/bin/bash

N="$1"
k="$2"
T="$3"

for i in {1..15}
do
  # nohup julia turn_flock_model.jl $N $k $T $i &
  nohup julia ../test/test_model.jl $N $k $T $i &
done
