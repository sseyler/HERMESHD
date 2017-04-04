#!/bin/bash

# Load OpenMPI and Intel Fortran compiler modules
module load openmpi/1.6.3/intel13.0/64


# Exports to set up paths Intel MKL

# Intel 17.0
# module load intel/17.0/fortran/64
# export MKLROOT=/nfs/packages/opt/Linux_x86_64/intel/17.0/mkl
# export MKL_TARGET_ARCH=intel64_lin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/${MKL_TARGET_ARCH}

# source $MKLROOT/bin/mklvars.sh intel64_lin

# Intel 13.0
export MKLROOT=/nfs/packages/opt/Linux_x86_64/intel/13.0/mkl
export MKL_TARGET_ARCH=intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/${MKL_TARGET_ARCH}

source $MKLROOT/bin/mklvars.sh intel64
