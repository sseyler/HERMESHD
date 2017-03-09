#!/usr/bin/make -f
##----------------------------------------------------------------------------
# Created with genmake.pl v1.1 on Sat Jan  5 21:18:00 2008
# genmake.pl home: http://muquit.com/muquit/software/
# Copryright: GNU GPL (http://www.gnu.org/copyleft/gpl.html)
##----------------------------------------------------------------------------

rm = /bin/rm -f
FORT = ifort
FORT = mpif90

DEFS =
PROGNAME = perseus_$(NAME)
TESTDIR = /scratch/sseyler/WORK/PERSEUSHD_testing
INCLUDES = -I. #-I$(MKLINCLUDE)
LIBS = lib_vtk_io.o mkl_vsl.o #-L$(MKLPATH)

DEFINES = $(INCLUDES) $(DEFS)
# FORTFLAGS = -O3 -coarray=distributed
# FORTFLAGS = -O3 -coarray -coarray-num-images=9
# FORTFLAGS = -O3 -coarray -coarray-num-images=4
# FORTFLAGS = -O3 -openmp
# FORTFLAGS = -O2 -xHost #-CB -fpe:0 -traceback -check noarg_temp_created #-fpp

##################################################################
# Intel compiler and MKL files
#-----------------------------------------------------------------
# ifortvars.sh -arch intel64 -platform linux
# IGNORE := $(shell bash -c "source ifortvars.sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
# include makeenv
MKLROOT = /nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/mkl
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64_lin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(MKLROOT)/lib/$(MKL_TARGET_ARCH)
# MKLBIN=$(MKLPATH)/bin
# MKLINCLUDE=$(MKLPATH)/include
MKL_TARGET_ARCH = intel64_lin
MISCFLAGS = -diag-disable 13003 # disable the warning for message verification

FORTFLAGS = -O2 -xHost -mkl -L$(MKLROOT)/lib/$(MKL_TARGET_ARCH) $(MISCFLAGS)
# FORTFLAGS = -O3 -xHost -mkl -L$(MKLROOT)/lib/$(MKL_TARGET_ARCH)
# FORTFLAGS = -O2 -xHost -CB -fpe:0 -traceback -check noarg_temp_created #-fpp
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
    OBJS = PERSEUSHydro_bitri_FH_test.f90
    MODULES =
endif

#********************************************
# Reference modular code
#********************************************
ifeq ($(NAME),mod)
    OBJS = perseushd_main.f90
    MODULES = parameters_mod.f90 ic_mod.f90 \
               prepare_time_advance_mod.f90 basis_funcs_mod.f90 \
               io_mod.f90
endif

#********************************************
# Semi-old modular code
#********************************************
ifeq ($(NAME),mod2)
    OBJS = main.f90
    MODULES = input.f90 parameters.f90 helpers.f90 \
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
    OBJS = main.f90
    MODULES =  input.f90 parameters.f90 basis_funcs.f90 helpers.f90 \
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
    OBJS = main.f90
    MODULES =  input_sod.f90 parameters_sod.f90 \
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
# 2D Poiseuille flow for development code
#********************************************
ifeq ($(NAME),poi2d)
    OBJS = main.f90
    MODULES =  input_2dpoiseuille.f90 parameters_2dpoiseuille.f90
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
	$(FORT) $(FORTFLAGS) -c $*.f90

all: $(PROGNAME)

$(PROGNAME) : $(OBJS)
	cp $(OBJS) $(MODULES) $(OTHERFILES) $(TESTDIR) && cd $(TESTDIR) && \
	$(FORT) $(FORTFLAGS) $(LIBS) $(MODULES) $(OBJS) -o $(PROGNAME)
#	strip $(PROGNAME)
#	$(PROGNAME)
	cd $(TESTDIR) && mpirun -n 16 $(PROGNAME) 2>&1 | tee $(PROGNAME).$(NAME).err

clean:
#	$(rm) $(PROGNAME) core *~
	$(rm) $(MODFILES)
##################################################################
