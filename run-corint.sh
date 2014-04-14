#!/bin/bash

cd /home/baxterma/run

for E in 1 2 5 7 10 15 20 25 30 40 50 60 70 80 100 200 500 1000 2000
do
   for z in 28 29 30 31 32 33 34
   do

      cat << eof > ./Ic-E$E-z$z.sh
#!/bin/bash
#PBS -N Ic-z$z-E$E
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=12
#PBS -r n
#PBS -V
#PBS -o ./runout/ic-z$z-E$E.out
#PBS -e ./runout/ic-z$z-E$E.err

cd /home/baxterma/run

./main $E $z > ./runout/res-z$z-E$E.txt
eof

      qsub Ic-E$E-z$z.sh

   done
done
