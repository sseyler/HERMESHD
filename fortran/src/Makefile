##----------------------------------------------------------------------------
# Created with genmake.pl v1.1 on Sat Jan  5 21:18:00 2008
# genmake.pl home: http://muquit.com/muquit/software/
# Copryright: GNU GPL (http://www.gnu.org/copyleft/gpl.html)
##----------------------------------------------------------------------------
rm=/bin/rm -f
FORT= ifort
FORT= mpif90

DEFS=
PROGNAME= perseus_$(NAME)
TESTDIR= /scratch/sseyler/WORK/PERSEUSHD_testing
INCLUDES= -I. #-I$(MKLINCLUDE)
LIBS= lib_vtk_io.o mkl_vsl.o #-L$(MKLPATH)

DEFINES= $(INCLUDES) $(DEFS)
# FORTFLAGS= -O3 -coarray=distributed
# FORTFLAGS= -O3 -coarray -coarray-num-images=9
# FORTFLAGS= -O3 -coarray -coarray-num-images=4
# FORTFLAGS= -O3 -openmp
# FORTFLAGS= -O2 -xHost #-CB -fpe:0 -traceback -check noarg_temp_created #-fpp

#################################
# ifortvars.sh -arch intel64 -platform linux
# IGNORE := $(shell bash -c "source ifortvars.sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
# include makeenv
MKLROOT=/nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/mkl
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64_lin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(MKLROOT)/lib/$(MKL_TARGET_ARCH)
# MKLBIN=$(MKLPATH)/bin
# MKLINCLUDE=$(MKLPATH)/include
MKL_TARGET_ARCH= intel64_lin
FORTFLAGS= -O2 -xHost -mkl -L$(MKLROOT)/lib/$(MKL_TARGET_ARCH)#-CB -fpe:0 -traceback -check noarg_temp_created #-fpp

# OBJS = dg_3D_hydro_bitri_H2O.f90
# OBJS = dg_3D_hydro_bitri_H2O_old.f90
# OBJS = dg_3D_hydro_bitri_version2_SEAN.f90
# OBJS = dg_3D_hydro_bitri_viscosity.f90

ifeq ($(NAME),bitri)
	OBJS = PERSEUSHydro_bitri_FH_test.f90
	MODS =
endif
ifeq ($(NAME),mod)
	OBJS = perseushd_main.f90
	MODS = parameters_mod.f90 initial_conditions_mod.f90 prepare_time_advance_mod.f90 auxiliary_mod.f90 basis_functions_mod.f90 io_mod.f90
endif

.c.o:
	$(rm) $@
	$(FORT) $(FORTFLAGS) -c $*.f90

all: $(PROGNAME)

$(PROGNAME) : $(OBJS)
	$(FORT) $(FORTFLAGS) $(LIBS) $(MODS) $(OBJS) -o $(PROGNAME)
#	strip $(PROGNAME)
#	$(PROGNAME)
	mv $(PROGNAME) $(TESTDIR) && cd $(TESTDIR) && mpiexec -n 16 $(PROGNAME)

clean:
	$(rm) $(OBJS) $(PROGNAME) core *~
