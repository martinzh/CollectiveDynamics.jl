#!/bin/bash

N="$1"
T="$2"

bash reps_vsk_2D.sh $N 0.001 $T
wait
bash reps_vsk_2D.sh $N 0.005 $T
wait
bash reps_vsk_2D.sh $N 0.0075 $T
wait
bash reps_vsk_2D.sh $N 0.01 $T
wait
bash reps_vsk_2D.sh $N 0.05 $T
wait
bash reps_vsk_2D.sh $N 0.075 $T
wait
bash reps_vsk_2D.sh $N 0.1 $T
wait
bash reps_vsk_2D.sh $N 0.3 $T
wait
bash reps_vsk_2D.sh $N 0.5 $T
