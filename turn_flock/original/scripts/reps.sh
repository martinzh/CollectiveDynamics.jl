#!/bin/bash

N="$1"
e="$2"
n="$3"
T="$4"

for i in {1..10}
do
  # nohup julia turn_flock_model.jl $N $k $T $i &
 time  nohup julia ../test/test_model.jl $N $e $n $T $i &
done
