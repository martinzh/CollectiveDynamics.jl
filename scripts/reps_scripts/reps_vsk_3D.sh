#!/bin/bash

N="$1"
k="$2"
T="$3"

for i in {1..10}
do
    time  nohup julia vicsek_3D.jl $N $k $T $i &
done
