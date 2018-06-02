#!/bin/bash

# Load cmake, OpenMPI, and Intel Fortran compiler modules
module load cmake/3.11/64
module load openmpi/1.6.3/intel13.0/64

export CC=mpicc
export CXX=mpicxx
export FC=mpif90

################################################################################
# Intel 13.0
export MKLROOT=/nfs/packages/opt/Linux_x86_64/intel/13.0/mkl
export MKL_TARGET_ARCH=intel64
export LD_LIBRARY_PATH=${MKLROOT}/lib/${MKL_TARGET_ARCH}:$LD_LIBRARY_PATH

source $MKLROOT/bin/mklvars.sh intel64
#-------------------------------------------------------------------------------
