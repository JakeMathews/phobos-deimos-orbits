#! /bin/bash

rm *.dat
rm *.out
gfortran simulation.f90
./a.out
xmgrace phobos-energy.dat