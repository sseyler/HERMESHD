#!/usr/bin/make -f
##----------------------------------------------------------------------------
# Created with genmake.pl v1.1 on Sat Jan  5 21:18:00 2008
# genmake.pl home: http://muquit.com/muquit/software/
# Copryright: GNU GPL (http://www.gnu.org/copyleft/gpl.html)
##----------------------------------------------------------------------------

rm = /bin/rm -f
FORT = ifort
MPIFORT = mpif90

DEFS =
PROGNAME = perseus_$(NAME)
TESTDIR = /scratch/sseyler/WORK/PERSEUSHD_testing

# NOTE: Seems to work when doing "include 'mpif.h'" in the code
MPIVERS = openmpi/1.6.3/intel13.0
MPIVERS = openmpi/2.1.0/intel17.0
INC_DIR = /nfs/packages/opt/Linux_x86_64/$(MPIVERS)/include

INCLUDES = -I. -I$(INC_DIR) #-I$(MKLINCLUDE)
OBJS = lib_vtk_io.o mkl_vsl.o #-L$(MKLPATH)

DEFINES = $(INCLUDES) $(DEFS)

##################################################################
# Intel compiler and MKL files
#-----------------------------------------------------------------
# ifortvars.sh -arch intel64 -platform linux
# IGNORE := $(shell bash -c "source ifortvars.sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
# include makeenv

# NOTE: libmkl_intel_lp64.so lives here: $(MKLROOT)/lib/intel64_lin
# NOTE: "mpif.h" is in /nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include

MKLROOT=/nfs/packages/opt/Linux_x86_64/intel/13.0/mkl
MKLPATH=$(MKLROOT)/lib
MKLINCLUDE=$(MKLROOT)/include
DYLD_LIBRARY_PATH=$(MKLPATH)

# MKLROOT = /nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/mkl
# MKL_TARGET_ARCH = intel64_lin
MISCFLAGS = -diag-disable 13003 # disable the warning for message verification

FORTFLAGS = -O2 -xHost -mkl -L$(MKLROOT)/lib/$(MKL_TARGET_ARCH) $(MISCFLAGS)
# FORTFLAGS = -O2 -xHost -CB -fpe:0 -traceback -check noarg_temp_created #-fpp
# FORTFLAGS = -O3 -coarray=distributed
# FORTFLAGS = -O3 -coarray -coarray-num-images=9
# FORTFLAGS = -O3 -coarray -coarray-num-images=4
# FORTFLAGS = -O3 -openmp
# FORTFLAGS = -O2 -xHost #-CB -fpe:0 -traceback -check noarg_temp_created #-fpp
##################################################################


##################################################################
# Source files
#-----------------------------------------------------------------
# OBJS = dg_3D_hydro_bitri_H2O.f90
# OBJS = dg_3D_hydro_bitri_H2O_old.f90
# OBJS = dg_3D_hydro_bitri_version2_SEAN.f90
# OBJS = dg_3D_hydro_bitri_viscosity.f90

OTHERFILES = mkl_vsl.mod mkl_vsl_type.mod mkl_vsl.o lib_vtk_io.mod lib_vtk_io.o

#********************************************
# Reference original code
#********************************************
ifeq ($(NAME),bitri)
    SRC = PERSEUSHydro_bitri_FH_test.f90
    MODSRC =
endif

#********************************************
# Reference modular code
#********************************************
ifeq ($(NAME),mod)
    SRC = perseushd_main.f90
    MODSRC = parameters_mod.f90 ic_mod.f90 \
               prepare_time_advance_mod.f90 basis_funcs_mod.f90 \
               io_mod.f90
endif

#********************************************
# Semi-old modular code
#********************************************
ifeq ($(NAME),mod2)
    SRC = main.f90
    MODSRC = input.f90 parameters.f90 helpers.f90 \
               init_conditions.f90 initialize.f90 \
               prep_time_advance.f90 basis_funcs.f90 io.f90
    MODFILES = input.mod parameters.mod helpers.mod \
               init_conditions.mod initialize.mod \
               prep_time_advance.mod basis_funcs.mod io.mod
endif

#********************************************
# Default modular (development) code
#********************************************
ifeq ($(NAME),modbc)
    SRC = main.f90
    MODSRC =  input.f90 parameters.f90 basis_funcs.f90 helpers.f90 \
               boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = input.mod parameters.mod basis_funcs.mod helpers.mod \
               boundary_defs.mod boundary_custom.mod boundary.mod \
               random.mod initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
ifeq ($(NAME),sod)
    SRC = main.f90
    MODSRC =  input_sod.f90 parameters_sod.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = input.mod parameters.mod \
               basis_funcs.mod helpers.mod \
               boundary_defs.mod boundary_custom.mod boundary.mod \
               random.mod initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
ifeq ($(NAME),sodpy)
    SRC = hermeshd.f90
    MODSRC =  LIB_VTK_IO.f90 input_sod.f90 parameters_sod.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = LIB_VTK_IO.mod input.mod parameters.mod \
               basis_funcs.mod helpers.mod \
               random.mod boundary_defs.mod boundary_custom.mod boundary.mod \
               initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
ifeq ($(NAME),sodpymkl)
    SRC = hermeshd.f90
    MODSRC =  LIB_VTK_IO.f90 mkl_vsl.f90 input_sod.f90 parameters_sod.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = LIB_VTK_IO.mod mkl_vsl_type.mod mkl_vsl.mod input.mod parameters.mod \
               basis_funcs.mod helpers.mod \
               random.mod boundary_defs.mod boundary_custom.mod boundary.mod \
               initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif


#********************************************
# 2D Poiseuille flow for development code
#********************************************
ifeq ($(NAME),poi2d)
    SRC = main.f90
    MODSRC =  input_2dpoiseuille.f90 parameters_2dpoiseuille.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = input.mod parameters.mod \
               basis_funcs.mod helpers.mod \
               boundary_defs.mod boundary_custom.mod boundary.mod \
               random.mod initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif
##################################################################


##################################################################
# Recipes
#-----------------------------------------------------------------
.c.o:
	$(rm) $@
	$(MPIFORT) $(FORTFLAGS) -c $*.f90

all: $(PROGNAME)

shlib:
	$(MPIFORT) -shared -fPIC $(FORTFLAGS) $(SRC) LIB_VTK_IO.f90 $(MODSRC) mkl_vsl.o -o $(PROGNAME).so

# NOTE: The compiled *.mod file for LIB_VTK_IO.f90 is lowercase and the object file is uppercase!!!
# NOTE: LIB_VTK_IO.f90 needs to be compiled BEFORE output.f90!!!

#env OMPI_FC=gfortran $(MPIFORT) -c -fPIC $(MODSRC) lib_vtk_io.mod  # LIB_VTK_IO.o lib_vtk_io.o mkl_vsl.o LIB_VTK_IO.f90
f2py:
	$(MPIFORT) -fc=$(FORT) -c -fPIC $(MODSRC) $(SRC)

$(PROGNAME) : $(SRC) $(MODSRC)
	cp $(SRC) $(MODSRC) $(OTHERFILES) $(TESTDIR) && cd $(TESTDIR) && \
	$(MPIFORT) $(FORTFLAGS) $(SRC) $(MODSRC) $(OBJS) -o $(PROGNAME)
	cd $(TESTDIR) && mpirun -n 16 $(PROGNAME) 2>&1 | tee $(PROGNAME).$(NAME).err

clean:
#	$(rm) $(PROGNAME) core *~
	$(rm) $(MODFILES) *.o .f2py_f2cmap f90wrap_*.f90 hermeshd.py* *.so
##################################################################
