read -p "File Name: " file1
filec=""$file1".c"
fileo=""$file1".o"
read -p ".Out File Name: " file2
fileout=""$file2".out"
read -p "Parameters Input: " params

gcc -Wall -I ~/gsl/include -c $filec 
gcc -L ~/gsl/lib $fileo -lgsl -lgslcblas -lm -o $fileout
./$fileout $params &


# 0.05 100 5E-4 1.25E-3 100 10000 1000 2400 20 922