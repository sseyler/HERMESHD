################################################################################
# CMake project file for hermeshd
################################################################################

cmake_minimum_required (VERSION 3.6.2)
project (test)
enable_language (Fortran)


list (APPEND CMAKE_MODULE_PATH cmake)

# this is where we will place the Fortran module files
set (SCRATCH /scratch/sseyler/WORK/HERMESHD_testing/)
set (EXECUTABLE_OUTPUT_PATH ${SCRATCH}/bin)
# set (CMAKE_Fortran_MODULE_DIRECTORY ${SCRATCH}/modules)

#-------------------------------------------------------------------------------


################################################################################
# make sure that the default is a RELEASE
#********************************************
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "None")
    # set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    #     "Choose the type of build, options are: None Debug Release."
    #     FORCE)
endif()
#-------------------------------------------------------------------------------


################################################################################
# Include CMake files
#********************************************
# these are paths that cmake will search for module files
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} cmake")
include(ConfigurePythonWrapping)
#-------------------------------------------------------------------------------


################################################################################
# Set compilation defaults
#********************************************
# Override defaults from the command line with:
#   cmake -DMAX_BUFFER=2048
#********************************************
# set(MAX_BUFFER "1024" CACHE STRING "Max buffer  make sure that the default is a RELEASEsize")
#
# # this is passed on to the code we compile
# add_definitions(-DMAX_BUFFER=${MAX_BUFFER})
#
# message(STATUS "Set max buffer to ${MAX_BUFFER}")
#-------------------------------------------------------------------------------


################################################################################
# MPI: default is ON
#********************************************
# Select this option on the command line with:
#   cmake -DENABLE_MPI=ON
#********************************************
# option(ENABLE_MPI "Enable MPI parallelization" ON)
#
# if (ENABLE_MPI)
#     message("MPI is enabled")
# else()
#     message("MPI is disabled")
# endif()

option(ENABLE_MKL "Enable Intel MKL (currently required)" ON)
if (ENABLE_MKL)
    message("Intel MKL is enabled")
    set(MKLROOT /nfs/packages/opt/Linux_x86_64/intel/13.0/mkl)
    set(MKL_TARGET_ARCH intel64)
    set(MKLPATH ${MKLROOT}/lib/${MKL_TARGET_ARCH})
else()
    message(WARNING "Intel MKL is disabled")
endif()
#-------------------------------------------------------------------------------


################################################################################
# default installation
#********************************************
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the build directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
#-------------------------------------------------------------------------------


################################################################################
# Set compiler-specific flags
#********************************************
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
    set(flags "-O2 -xHost -mkl -diag-disable 13003")
endif()
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} ${flags}")
set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${flags}")
#-------------------------------------------------------------------------------


################################################################################
# Compile
#********************************************
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB sources ${SOURCE_DIR}/*.f90)

# Sources for normal compilation of executable
# NOTE: It will be useful to have separate subdirectories, probably
set(default_sources ${sources})
list(REMOVE_ITEM default_sources ${SOURCE_DIR}/hermeshd.f90)
list(REMOVE_ITEM default_sources ${SOURCE_DIR}/f90wrap_hermeshd.f90)
list(REMOVE_ITEM default_sources ${SOURCE_DIR}/spatial.f90)

set (PROG "hermeshd")

#-------------------------------------------------------------------------------


option(ENABLE_PYWRAP "Enable Python wrapping and module compilation" OFF)
if (ENABLE_PYWRAP)
    message("Python module compilation is enabled")

    set (PYWRAP_NAME ${PROG})

    # Sources for Python module compilation
    set(module_sources ${sources})
    list(REMOVE_ITEM module_sources ${SOURCE_DIR}/main.f90)
    list(REMOVE_ITEM module_sources ${SOURCE_DIR}/spatial.f90)
    list(REMOVE_ITEM module_sources ${SOURCE_DIR}/f90wrap_hermeshd.f90)

    # Compile sources into static library (with -fPIC flag)
    # Used by _${PYWRAP_MOD}.so to generate the Python extension module
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${flags} -fPIC")
    add_library(${PYWRAP_NAME} ${module_sources})

    # Generates new target _${PYWRAP_NAME} that depends on ${PYWRAP_NAME}
    PYWRAP(
        MOD ${PYWRAP_NAME}
        SRC ${module_sources}
        EXT_LIBS "-L${MKLPATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx"
    )

    # Compilation depends on
    add_dependencies(_${PYWRAP_NAME} ${PYWRAP_NAME})

    # Generate executable from same static library as used by PYWRAP?
    # add_executable(${PROG} main.f90)
    # target_link_libraries(${PROG} ${PYWRAP_NAME})

else()
    message("Normal compilation of executable will take place")

    add_executable(${PROG} default_sources)
endif()

# install (TARGETS _${F2PY_ADD_MODULE}.so RUNTIME DESTINATION ${SCRATCH})
# install (TARGETS ${PYFILES} RUNTIME DESTINATION ${SCRATCH})
