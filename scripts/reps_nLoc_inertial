#!/bin/bash

N="$1"
e="$2"
T="$3"
n="$4"
t="$5"

for i in {1..10}
do
    time  nohup julia inertial_nonLocal_simScript.jl $N $e $T $n $t $i &
done
