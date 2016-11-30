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
source /nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/mkl/bin/mklvars.sh intel64 mod
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(MKLROOT)/lib/$(MKL_TARGET_ARCH)
```

Where ``$(MKLROOT)`` is set by the *mklvars.sh* script and ``$(MKL_TARGET_ARCH)``
should be set manually to "*intel64_lin*" (or "*intel64*", which is linked to
"*intel64_lin*").

## Try compiling manually

Note that the `use` statement in a fortran file requires that the corresponding
*.mod* file be included with `-I/path/to/dir/with/mod/file` when compiling. To
use a module in *my_prog.f90* called `module_name` (i.e. `use MODULE_NAME`) in
a file *my_module.f90*, first compile the module with something like

```bash
ifort -c my_module.f90
```

This will produce a module file with the module's name *MODULE_NAME.mod* and an
object file *my_module.o*. Then, compile the program with

```bash
ifort -c my_prog.f90 -I/path/to/dir/with/mod/file
```

where the `-I<dir>` options provides the compiler with the path to search and
locate the *MODULE_NAME.mod* file.



### MPI compilation

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
