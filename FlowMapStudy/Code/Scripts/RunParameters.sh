#!/bin/bash

Nx=3
Ny=3
Nz=2
PRINT_VTK=0
DO_FTLE=1
ADVECT_PARTICLES=7
Samples_cubic=10

for ENDTIME in 1 2 4 8 16 32 64 128
do
	for STEPSIZE in 0.01 0.001 0.0001 0.00001
	do
		for AVG in 0 1 2 3
		do
			echo $Nx $Ny $Nz $ENDTIME $PRINT_VTK $DO_FTLE $ADVECT_PARTICLES $STEPSIZE $Samples_cubic $AVG
			../bin/runTokamak_small.out $Nx $Ny $Nz $ENDTIME $PRINT_VTK $DO_FTLE $ADVECT_PARTICLES $STEPSIZE $Samples_cubic &> ../bin/OUT\_SCRIPT/Tokamak\_Small\_${Nx}\_${Ny}\_${Nz}\_${ENDTIME}\_${PRINT_VTK}\_${DO_FTLE}\_${ADVECT_PARTICLES}\_${STEPSIZE}\_${Samples_cubic}\_${AVG}
			../bin/runTokamak.out $Nx $Ny $Nz $ENDTIME $PRINT_VTK $DO_FTLE $ADVECT_PARTICLES $STEPSIZE $Samples_cubic &> ../bin/OUT\_SCRIPT/Tokamak\_Regular\_${Nx}\_${Ny}\_${Nz}\_${ENDTIME}\_${PRINT_VTK}\_${DO_FTLE}\_${ADVECT_PARTICLES}\_${STEPSIZE}\_${Samples_cubic}\_${AVG}
		done
	done
done



