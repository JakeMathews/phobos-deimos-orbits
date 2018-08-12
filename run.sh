#! /bin/bash

rm *.dat
rm *.out
gfortran simulation.f90
./a.out
xmgrace deimos.dat