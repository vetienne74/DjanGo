#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take a few minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

./clean_all

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/
cp ../Par_file .
cp ../interfaces_elastic_analytic.dat .
cp ../SOURCE SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# runs database generation
echo
echo "  running mesher..."
echo

export OMP_NUM_THREADS=2

# GNU
export GOMP_CPU_AFFINITY=0,1
mpirun -n 1 ${SPECFEM2D_DIR}/bin/xmeshfem2D

# runs simulation
echo
echo "  running solver..."
echo
mpirun -n 1 ${SPECFEM2D_DIR}/bin/xspecfem2D

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

cp OUTPUT_FILES/Ux_file_single.bin ./vx.time.rec.specfem2d.out
cp OUTPUT_FILES/Uz_file_single.bin ./vz.time.rec.specfem2d.out

#./clean_all
exit

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
