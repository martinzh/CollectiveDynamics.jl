#!/bin/bash

N="$1"
T="$2"

bash reps_locNloc_top_2D.sh $N 0.0 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 0.01 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 0.025 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 0.05 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 0.1 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 0.25 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 0.5 0.5 $T
wait
bash reps_locNloc_top_2D.sh $N 1.0 0.5 $T
wait
