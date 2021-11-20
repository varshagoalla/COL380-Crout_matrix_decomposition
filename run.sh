#!/bin/bash

input_size=$1
input_file=$2
num_threads=$3
strategy=$4

if [ $4 == 0 ] 
then
	./g $1 $2 $3 $4
elif [ $4 == 1 ] 
then
	./g $1 $2 $3 $4
elif [ $4 == 2 ] 
then
	./g $1 $2 $3 $4
elif [ $4 == 3 ] 
then
	./g $1 $2 $3 $4
elif [ $4 == 4 ] 
then
	mpiexec -n $3 ./m $1 $2 $3 $4
else
	echo "wrong input strategy"
fi
