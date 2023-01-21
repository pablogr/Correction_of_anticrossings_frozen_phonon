#!/bin/bash

rm FPTripleAnticrossing.x

#gfortran-mp-7 prueba.f90 -L/usr/local/lib/lapack -llapack -o prueba.x
gfortran-mp-7 FPTripleAnticrossing.f90  -llapack -o FPTripleAnticrossing.x

./FPTripleAnticrossing.x
