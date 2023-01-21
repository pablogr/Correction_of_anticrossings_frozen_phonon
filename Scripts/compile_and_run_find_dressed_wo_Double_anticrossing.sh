#!/bin/bash

rm FPDoubleAnticrossing.x

gfortran-mp-7 FPDoubleAnticrossing.f90 -llapack -o FPDoubleAnticrossing.x

./FPDoubleAnticrossing.x
