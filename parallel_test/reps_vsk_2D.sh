#!/bin/bash

N="$1"
r="$2"
T="$3"

for i in {1..5}
do
    time nohup julia vicsek_2D.jl $N $r $T $i &
done
