#!/bin/bash

N="$1"
T="$2"

bash reps.sh $N 0 $T
bash reps.sh $N 2 $T
bash reps.sh $N 5 $T
bash reps.sh $N 7 $T
