## Load modules:

### Load the Intel Fortran compiler

```bash
module load intel/17.0/fortran/64
```

### This seems to be needed for licensing...

```bash
module load intel/13.0/64
```

### Load an MPI library

OpenMPI works, but try Intel's MPI library at some point

```bash
openmpi/1.6.3/intel13.0/64
```


## Load environment variables

```bash
/nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/mkl/bin/mklvars.sh intel64 mod
```

## Try compiling manually with

```bash
mpif90 -O2 -xHost lib_vtk_io.o dg_3D_hydro_bitri_version2_FH_test.f90 -o perseus_fh_test
```

```bash
mpif90 -O2 -xHost -mkl -I. -I${MKLROOT}/include lib_vtk_io.o dg_3D_hydro_bitri_version2_FH_test.f90 -o perseus_fh_test
```

## Run with either

```bash
mpiexec -n 16 perseus_fh_test
```

```bash
mpirun -n 16 perseus_fh_test
```


# Compiling with Intel MKL

Might need these at top of PERSEUS code:

```fortran
! include '/nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/mkl/include/mkl_vsl.f90'
include 'mkl_vsl.f90'

! Get Intel MKL types and functions
use MKL_VSL_TYPE
use MKL_VSL
```
