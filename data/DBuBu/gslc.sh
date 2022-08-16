# read -p "File Name: " file1
# filec=""$file1".c"
# fileo=""$file1".o"
# read -p ".Out File Name: " file2
# fileout=""$file2".out"
# read -p "Parameters Input: " params
file1="hill"
filec=""$file1".c"
fileo=""$file1".o"
file2="t"
fileout=""$file2".out"
params="0.05 100 5E-4 1.25E-3 100 10000 1000 2400 490 "

gcc -Wall -I ~/gsl/include -c $filec 
gcc -L/usr/local/Cellar/gsl/2.5/lib $fileo -lgsl -lgslcblas -lm -o $fileout
./$fileout $params  

# 0.05 100 5E-4 2.5E-1 1 100 1000 300 409
