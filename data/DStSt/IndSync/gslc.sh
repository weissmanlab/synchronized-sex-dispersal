file1="hill_is"
filec=""$file1".c"
fileo=""$file1".o"
file2="t"
fileout=""$file2".out"

gcc -Wall -I ~/gsl/include -c $filec 
gcc -L/usr/local/Cellar/gsl/2.5/lib $fileo -lgsl -lgslcblas -lm -o $fileout
params="0.05 100 5E-4 1.25E-3 100 10 1000 2400 922"
./$fileout $params  &
params="0.05 100 5E-4 1.25E-3 100 100 1000 2400 922"
./$fileout $params  &
params="0.05 100 5E-4 1.25E-3 100 1000 1000 2400 922"
./$fileout $params  &
params="0.05 100 5E-4 1.25E-3 100 10000 1000 2400 922"
./$fileout $params  &

# 0.05 1 5E-4 1.25E-3 100 10000 1000 2400 924