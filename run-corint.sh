#!/bin/bash

# $1 = name of calculation

mkdir -p runout

echo -e "\nStart: $(date)\n"
start=$(date +%s)

for z in 28 29 30 31 32 33 34
do

   mkdir -p ./output/$1-un/z$z
   mkdir -p ./output/$1/z$z

   for E in 1 2 5 7 10 15 20 25 30 40 50 60 70 80 100 200 500 1000 2000
   do
   
      sed -i "15s/.*/   int E = $E;/" main.cpp
      sed -i "16s/.*/   int z = $z;/" main.cpp
      sed -i "28s/.*/   ofstream writefile  ( "./output/$1-un/z$z/$1-E$E-z$z.txt" );/" main.cpp

      make -sj

      ./main > ./runout/E$E-z$z-times.txt

      echo "done E: $E  z: $z"
      
      sed -i "59s/.*/\treadpath  = ".output/$1-un/z$z/$1-E$E-z$z.txt"/" columns.py
      sed -i "60s/.*/\treadpath  = ".output/$1/z$z/$1-E$E-z$z.txt"/" columns.py
      
      python columns.py

   done
done

echo -e "End: $(date)\n"
end=$(date +%s)

echo "Done in $(($end-$start)) s"

sed -i "38s/.*/folder  = "$1"/" average.py

python average.py

rm -r ./output/$1-un/z$z
