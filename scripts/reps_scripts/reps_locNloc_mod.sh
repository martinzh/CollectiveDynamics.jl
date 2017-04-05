#!/bin/bash

N="$1"
k="$2"
w="$3"
T="$4"
d="$5"
m="$6"

for i in {1..10}
do
    time  nohup julia locNloc_mod_simScript.jl $N $k $w $T $d $m $i &
done
