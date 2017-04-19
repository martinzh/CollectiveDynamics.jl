#!/bin/bash

N="$1"
k="$2"
w="$3"
T="$4"

for i in {1..10}
do
    time  nohup julia 3D_locNloc_simScript_top.jl $N $k $w $T $i &
done
