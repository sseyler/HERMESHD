#!/usr/bin/make -f
##----------------------------------------------------------------------------
# Created with genmake.pl v1.1 on Sat Jan  5 21:18:00 2008
# genmake.pl home: http://muquit.com/muquit/software/
# Copryright: GNU GPL (http://www.gnu.org/copyleft/gpl.html)
##----------------------------------------------------------------------------

rm = /bin/rm -f
F90 = ifort
MPIF90 = mpif90

PROGNAME = perseus_$(NAME)
RUNDIR = /scratch/sseyler/WORK/PERSEUSHD_testing

# NOTE: Seems to work when doing "include 'mpif.h'" in the code
# MPIVERS = openmpi/1.6.3/intel13.0
# MPIVERS = openmpi/2.1.0/intel17.0
# INC_DIR = /nfs/packages/opt/Linux_x86_64/$(MPIVERS)/include

INCLUDES = -I. -I$(INC_DIR)
OBJS = LIB_VTK_IO.o mkl_vsl.o

##################################################################
# Intel compiler and MKL files
#-----------------------------------------------------------------
# ifortvars.sh -arch intel64 -platform linux
# IGNORE := $(shell bash -c "source ifortvars.sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
# include makeenv

# NOTE: libmkl_intel_lp64.so lives here: $(MKLROOT)/lib/intel64_lin
# NOTE: "mpif.h" is in /nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include

# MKLROOT = /nfs/packages/opt/Linux_x86_64/intel/13.0/mkl
MKLROOT = /nfs/packages/opt/Linux_x86_64/intel/17.0/mkl
MKLPATH = $(MKLROOT)/lib/intel64
MKLINCLUDE = $(MKLROOT)/include
DYLD_LIBRARY_PATH = $(MKLPATH)

MKL_TARGET_ARCH = intel64_lin
MISCFLAGS = -diag-disable 13003 # disable the warning for message verification

# F90FLAGS = -O2 -xHost -mkl $(MISCFLAGS) -I$(INCLUDES) -I$(MKLINCLUDE) -L$(MKLPATH)/$(MKL_TARGET_ARCH)
F90FLAGS = -O2 -xHost -mkl $(MISCFLAGS)
##################################################################


ifndef $(NPROC)
	NPROC = 16
endif



##################################################################
# Source files
#-----------------------------------------------------------------

#********************************************
# Reference original code
#********************************************
ifeq ($(NAME),bitri)
    SRC = PERSEUSHydro_bitri_FH_test.f90
    MODSRC =
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
    MODSRC =  input.f90 params.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = input.mod params.mod \
               basis_funcs.mod helpers.mod \
               boundary_defs.mod boundary_custom.mod boundary.mod \
               random.mod initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
ifeq ($(NAME),sodmkl)
    SRC = main.f90
    MODSRC =  LIB_VTK_IO.f90 mkl_vsl.f90 input.f90 params.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = LIB_VTK_IO.mod mkl_vsl_type.mod mkl_vsl.mod input.mod params.mod \
               basis_funcs.mod helpers.mod \
               random.mod boundary_defs.mod boundary_custom.mod boundary.mod \
               initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
ifeq ($(NAME),sodpy)
    SRC = hermeshd.f90
    MODSRC =  LIB_VTK_IO.f90 input.f90 params.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = LIB_VTK_IO.mod input.mod params.mod \
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
    MODSRC =  LIB_VTK_IO.f90 mkl_vsl.f90 input.f90 params.f90 \
               basis_funcs.f90 helpers.f90 boundary.f90 \
               random.f90 initialcon.f90 initialize.f90 \
               prepare_step.f90 sources.f90 flux.f90 \
               integrator.f90 output.f90
    MODFILES = LIB_VTK_IO.mod mkl_vsl_type.mod mkl_vsl.mod input.mod params.mod \
               basis_funcs.mod helpers.mod \
               random.mod boundary_defs.mod boundary_custom.mod boundary.mod \
               initialcon.mod initialize.mod \
               prepare_step.mod sources.mod flux.mod \
               integrator.mod output.mod
endif
##################################################################


##################################################################
# Recipes
#-----------------------------------------------------------------
# .c.o:
# 	$(rm) $@
# 	$(MPIF90) $(F90FLAGS) -c $*.f90
#
# all: $(PROGNAME)
#
# shlib:
# 	$(MPIF90) -shared -fPIC $(F90FLAGS) $(SRC) LIB_VTK_IO.f90 $(MODSRC) mkl_vsl.o -o $(PROGNAME).so

# $(PROGNAME) : $(SRC) $(MODSRC)
# 	cp $(SRC) $(MODSRC) $(OBJS) $(RUNDIR) && cd $(RUNDIR) && \
# 	$(MPIF90) $(F90FLAGS) $(SRC) $(MODSRC) $(OBJS) -o $(PROGNAME)
# 	cd $(RUNDIR) && mpirun -n 16 $(PROGNAME) 2>&1 | tee $(PROGNAME).$(NAME).err

f2py:
	$(MPIF90) -fc=$(F90) -c -fPIC $(MODSRC) $(SRC)

run: $(RUNDIR)/$(PROGNAME)
	cd $(RUNDIR) && mpirun -n $(NPROC) $(PROGNAME) 2>&1 | tee $(PROGNAME).$(NAME).err

$(RUNDIR)/$(PROGNAME): $(MODSRC) $(SRC) cprun
	cd $(RUNDIR) && $(MPIF90) $(F90FLAGS) $(MODSRC) $(SRC) -o $(PROGNAME)

cprun:
	cp $(SRC) $(MODSRC) $(RUNDIR)

clean:
#	$(rm) $(PROGNAME) core *~
	$(rm) $(MODFILES) *.o .f2py_f2cmap f90wrap_*.f90 hermeshd.py* *.so
##################################################################
