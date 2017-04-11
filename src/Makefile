#!/usr/bin/make -f
rm = /bin/rm -f

################################################################################
# Definitions
#-----------------------------------------------------------------

#--- Fortran compilers -----------
F90 = ifort
MPIF90 = mpif90

ifndef $(NPROC)
	NPROC = 16
endif

ifndef $(BUILDDIR)
	BUILDDIR = build
endif

ifndef $(EXEC)
	EXEC = perseus
endif

ifndef $(RUNDIR)
	# RUNDIR = /scratch/sseyler/WORK/PERSEUSHD_testing
	RUNDIR = /scratch/sseyler/WORK/HERMESHD_testing
endif
#-----------------------------------------------------------------


# NOTE: Seems to work when doing "include 'mpif.h'" in the code
# MPIVERS = openmpi/1.6.3/intel13.0
# MPIVERS = openmpi/2.1.0/intel17.0
# INC_DIR = /nfs/packages/opt/Linux_x86_64/$(MPIVERS)/include

INCLUDES = -I. -I$(INC_DIR)

################################################################################
# Intel compiler and MKL files
#-----------------------------------------------------------------
# ifortvars.sh -arch intel64 -platform linux
# IGNORE := $(shell bash -c "source ifortvars.sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
# include makeenv

# NOTE: libmkl_intel_lp64.so lives here: $(MKLROOT)/lib/intel64_lin
# NOTE: "mpif.h" is in /nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include

MISCFLAGS = -diag-disable 13003 # disable the warning for ifort message verification
# F90FLAGS = -O2 -xHost -mkl $(MISCFLAGS) -I$(INCLUDES) -I$(MKLINCLUDE) -L$(MKLPATH)/$(MKL_TARGET_ARCH)
F90FLAGS = -O2 -xHost -mkl $(MISCFLAGS)
################################################################################


################################################################################
# Source files
#-----------------------------------------------------------------

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
SRC      = main.f90
MODSRC   =  LIB_VTK_IO.f90 mkl_vsl.f90 \
			input.f90 params.f90 basis_funcs.f90 helpers.f90 random.f90 \
			boundary.f90 initialcon.f90 initialize.f90 prepare_step.f90 \
			sources.f90 flux.f90 integrator.f90 output.f90
MODFILES = LIB_VTK_IO.mod mkl_vsl_type.mod mkl_vsl.mod \
			input.mod params.mod basis_funcs.mod helpers.mod random.mod \
			boundary_defs.mod boundary_custom.mod boundary.mod \
			initialcon.mod initialize.mod prepare_step.mod sources.mod flux.mod \
			integrator.mod output.mod

#********************************************
# Sod Shock Tube 1D for development code
#********************************************
ifeq ($(VERS),pymkl)
	EXEC    = hermes
	NAME    = hermeshd
	SRC     = $(NAME).f90
	F90WRAPSRC = f90wrap_$(NAME).f90

	PYMAIN    = test_$(NAME).py
	PYSOD     = test_sod.py
	PYWRAPSRC = $(NAME).py
	PYSHARED  = _$(NAME).so
	PYOBJS    = $(PYMAIN) $(PYWRAPSRC) $(PYSHARED) $(PYSOD)
	# PYRUN = $(patsubst %, $(RUNDIR)/%, $(PYOBJS))

	OBJFILES = $(MODSRC:.f90=.o) $(NAME).o  # Get list of object files to be produced

	# Explicit MKL paths needed when using f2py
	MKLROOT = /nfs/packages/opt/Linux_x86_64/intel/13.0/mkl
	MKL_TARGET_ARCH = intel64
	# MKLROOT = /nfs/packages/opt/Linux_x86_64/intel/17.0/mkl
	# MKL_TARGET_ARCH = intel64_lin

	MKLPATH = $(MKLROOT)/lib/$(MKL_TARGET_ARCH)
	MKLINCLUDE = $(MKLROOT)/include
	DYLD_LIBRARY_PATH = $(MKLPATH)
endif

# OBJECTS = $(SOURCES:f90=o)
# OBJECTS = $(patsubst $(SOURCEDIR)/%.f90, $(BUILDDIR)/%.o, $(SOURCES))
################################################################################


################################################################################
# For printing out information
#-----------------------------------------------------------------
LINE1 = "--------------------------------------------------------------------"
LINE2 = "===================================================================="
LINE1N = \n$(LINE1)\n
LINE2N = \n$(LINE2)\n

STAGE1 = "\n$(LINE2N) Stage 1: generating Python wrappers for $(SRC)... $(LINE1N)"
STAGE2 = "\n$(LINE2N) Stage 2: compiling Fortran modules, generating object files... $(LINE1N)"
STAGE3 = "\n$(LINE2N) Stage 3: building extension module $(PYSHARED)... $(LINE1N)"


################################################################################
# Fortran recipes
#-----------------------------------------------------------------
all: $(BUILDDIR)/$(EXEC) run

run:
	@echo "$(LINE2N) Running $(EXEC)... $(LINE1N)"
	cp $(BUILDDIR)/$(EXEC) $(RUNDIR)
	cd $(RUNDIR) && mpirun -n $(NPROC) $(EXEC) 2>&1 | tee out.$(EXEC).run

$(BUILDDIR)/$(EXEC): $(MODSRC) $(SRC) cp-src
	cd $(BUILDDIR) && $(MPIF90) $(F90FLAGS) $(MODSRC) $(SRC) -o $(EXEC)

cp-src: $(MODSRC) $(SRC) | $(BUILDDIR)
	@echo "\n>>> Copying source files to build directory..."
	cp $(MODSRC) $(SRC) $(BUILDDIR)


#-----------------------------------------------------------------
# Python recipes
#-----------------------------------------------------------------
pyall: wrap pycompile pybuild

pyrun: cp-py-run
	@echo "$(LINE2N) Running $(PYMAIN)... $(LINE1N)"
	@echo "\n>>> Changing to run directory..."
	cd $(RUNDIR) && mpirun -n $(NPROC) python $(PYMAIN) 2>&1 | tee out.$(NAME).pyrun

pybuild: $(patsubst %, $(BUILDDIR)/%, $(OBJFILES)) $(BUILDDIR)/$(F90WRAPSRC)
	@echo $(STAGE3)
	cd $(BUILDDIR) && f2py-f90wrap --f90exec=$(MPIF90) --opt="-O2"				\
		--f90flags="-fc=ifort $(F90FLAGS)" -c -m _$(NAME) $(OBJFILES) $(F90WRAPSRC)				\
		-L$(MKLPATH) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx
	@echo "\nStage 3 completed. $(LINE1N)"

pybuild-lib: $(patsubst %, $(BUILDDIR)/%, $(OBJFILES)) $(BUILDDIR)/$(F90WRAPSRC)
	@echo $(STAGE3)
	cd $(BUILDDIR) && f2py-f90wrap --f90exec=$(MPIF90) --opt="-O2"				\
		--f90flags="-fc=ifort $(F90FLAGS)" -c -m _$(NAME) $(F90WRAPSRC)				\
		-L. -lhermeshd -L$(MKLPATH) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx
	@echo "\nStage 3 completed. $(LINE1N)"

pycompile: cp-py-bld cp-src
	@echo $(STAGE2)
	cd $(BUILDDIR) && $(MPIF90) $(F90FLAGS) -fc=$(F90) -c -fPIC $(MODSRC) $(SRC)
	@echo "\nStage 2 completed. $(LINE1N)"
	@echo $(OBJFILES)

pywrap:
	@echo $(STAGE1)
	f90wrap -m $(NAME) $(SRC) --only main setup step cleanup generate_output
	@echo "\nStage 1 completed. $(LINE1N)"

cp-py-bld: $(F90WRAPSRC) | $(BUILDDIR)
	@echo "\n>>> Copying Python files to build directory..."
	cp $(F90WRAPSRC) $(PYWRAPSRC) $(PYMAIN) $(PYSOD) $(BUILDDIR)

cp-py-run:
	@echo "\n>>> Copying Python files to run directory..."
	cp $(PYOBJS:%=$(BUILDDIR)/%) $(RUNDIR)

#-----------------------------------------------------------------
$(BUILDDIR):
	mkdir -p $@

clean:
	$(rm) -r $(BUILDDIR)

clean-src:
	$(rm) $(MODFILES) *.o .f2py_f2cmap $(WRAPSRC) $(PYWRAPSRC)* $(PYSHARED) lib_vtk_io.mod
################################################################################


##################################################################
# Recipes
#-----------------------------------------------------------------
# .c.o:
# 	$(rm) $@
# 	$(MPIF90) $(F90FLAGS) -c $*.f90
#
# all: $(EXEC)
#
# shlib:
# 	$(MPIF90) -shared -fPIC $(F90FLAGS) $(SRC) LIB_VTK_IO.f90 $(MODSRC) mkl_vsl.o -o $(EXEC).so

# $(EXEC) : $(SRC) $(MODSRC)
# 	cp $(SRC) $(MODSRC) $(OBJS) $(RUNDIR) && cd $(RUNDIR) && \
# 	$(MPIF90) $(F90FLAGS) $(MODSRC) $(SRC) $(OBJS) -o $(EXEC)
# 	cd $(RUNDIR) && mpirun -n 16 $(EXEC) 2>&1 | tee $(EXEC).$(VERS).out
