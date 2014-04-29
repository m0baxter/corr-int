#!/bin/bash

for E in 1 2 5 7 10 15 20 25 30 40 50 60 70 80 100 200 500 1000 2000
do
   for z in 28 29 30 31 32 33 34
   do
      rm Ic-E$E-z$z.sh
   done
done
