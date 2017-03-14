#!/bin/bash

N="$1"
e="$2"
T="$3"
t="$4"

for i in {1..10}
do
    time  nohup julia inertial_simScript.jl $N $e $T $t $i &
done
