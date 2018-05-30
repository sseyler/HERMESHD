################################################################################
# CMake project file for hermeshd
################################################################################

cmake_minimum_required (VERSION 3.6.2)
project (test)
enable_language (Fortran)

################################################################################
# default installation
#********************************************
# Location of CMake modules (*.cmake)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# get_filename_component (default_prefix ".." ABSOLUTE)
# set (default_prefix ${CMAKE_CURRENT_SOURCE_DIR})
# set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#       "Choose the installation directory; by default it installs in the current directory."
#       FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set (SCRATCH "/scratch/sseyler/WORK/BWAnnualReport2018/")
set (CMAKE_INSTALL_PREFIX ${SCRATCH})
# set (EXECUTABLE_OUTPUT_PATH ${SCRATCH}/bin)
# set (LIBRARY_OUTPUT_PATH ${SCRATCH}/lib)

# Location to place Fortran module files (*.mod)
set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
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
# if (ENABLE_MPI)
#     message("MPI is enabled")
# else()
#     message("MPI is disabled")
# endif()


option(MKL_PATH "Path to Intel MKL libraries" nfs/packages/opt/Linux_x86_64/intel/13.0/mkl/lib/intel64)

# option(MKL_LIBS "Intel MKL libraries" mkl_intel_lp64 mkl_sequential mkl_core mkl_vml_avx)
set (MKL_LIBS mkl_intel_lp64 mkl_sequential mkl_core mkl_vml_avx)
# set (MKL_LIBS "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx")

# option(ENABLE_MKL "Enable Intel MKL (currently required)" ON)
# if (ENABLE_MKL)
#     message("Intel MKL is enabled")
#     set(MKLROOT /nfs/packages/opt/Linux_x86_64/intel/13.0/mkl)
#     set(MKL_TARGET_ARCH intel64)
#     set(MKLPATH ${MKLROOT}/lib/${MKL_TARGET_ARCH})
# else()
#     message(WARNING "Intel MKL is disabled (currently required)")
# endif()
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
# list(REMOVE_ITEM default_sources ${SOURCE_DIR}/spatial.f90)

set (PROG "hermes")

# Sources for Python module compilation
set(module_sources ${sources})
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/main.f90)
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/f90wrap_hermeshd.f90)
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/reservoir.f90)
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/flux.f90)
# list(REMOVE_ITEM module_sources ${SOURCE_DIR}/spatial.f90)

set (LIB_NAME "hermeshd")
#-------------------------------------------------------------------------------

# Compile sources into static library (with -fPIC flag)
# Used by _${PYWRAP_MOD}.so to generate the Python extension module
message ("Will compile static library lib${PYWRAP_NAME}.a with the flags: ${CMAKE_Fortran_FLAGS}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")  # "-c" flag implied by "add_library"
add_library(${LIB_NAME} ${module_sources})

option(ENABLE_PYWRAP "Enable Python wrapping and module compilation" OFF)
if (ENABLE_PYWRAP)
    message("Python module compilation is enabled")

    set(PYMOD_NAME ${LIB_NAME})

    # Generates new target _${PYWRAP_NAME} that depends on ${PYWRAP_NAME}
    PYWRAP(
        MOD ${PYMOD_NAME}
        SRC ${module_sources}
        SUBPROGRAMS "main setup step cleanup generate_output"
        LIBS ${LIB_NAME}
        EXT_LIBS ${MKL_LIBS}
        EXT_LIBS_LOC ${MKL_PATH}
    )
    # install (FILES ${PROJECT_BINARY_DIR}/_${PYWRAP_NAME}.so ${PROJECT_BINARY_DIR}/${PYWRAP_NAME}.py DESTINATION lib)
else()
    message("Normal compilation of executable will take place")
    # install (TARGETS ${PROG} RUNTIME DESTINATION bin)
endif()

add_executable(${PROG} ${SOURCE_DIR}/main.f90)
target_link_libraries(${PROG} ${LIB_NAME})
