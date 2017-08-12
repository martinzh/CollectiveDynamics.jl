#!/bin/bash

N="$1"
k="$2"
w="$3"
T="$4"

for i in {1..10}
do
    time  nohup julia 2D_locNloc_simScript_met.jl $N $k $w $T $i &
done

wait
