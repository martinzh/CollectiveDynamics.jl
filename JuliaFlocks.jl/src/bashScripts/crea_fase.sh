#!/bin/bash

for i in {0..19}
do

  var=$(echo "0.005*$i" | bc -l)
  # echo cls_data0"$var".dat
  # echo awk '{print $2}' cls_data0$var.dat | tail -1 > tmp
  awk '{print $2}' cls_data0$var.dat | tail -1 > tmp
  echo "0$var" > tmp2
  paste tmp2 tmp >> fase.dat
done

  # echo awk '{print $2}' cls_data1.0.dat | tail -1 > tmp
  awk '{print $2}' cls_data1.0.dat | tail -1 > tmp
  echo 1.0 > tmp2
  paste tmp2 tmp >> fase.dat

  rm tmp tmp2
