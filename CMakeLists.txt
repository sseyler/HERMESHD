################################################################################
# CMake project file for hermeshd
################################################################################

cmake_minimum_required(VERSION 3.6.2)
project (test)
enable_language (Fortran)

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)


# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# Set compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
    set(flags "-O2 -xHost -mkl -diag-disable 13003")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} ${flags}")



#
# Compile.
#
file(GLOB sources src/*.f90)

set (EXECUTABLES "hermes")

add_executable(hermes ${sources})

install (TARGETS ${EXECUTABLES}
         RUNTIME DESTINATION "bin")




################################################################################
# Include any header files
# include_directories(include)

#C an manually add the sources using the set command as follows:
# set(SOURCES src/mainapp.cpp src/Student.cpp)

# However, the file(GLOB...) allows for wildcard additions:
# file(GLOB SOURCES "src/*.f90")
#
# add_executable(perseus ${SOURCES})
#-------------------------------------------------------------------------------



################################################################################
# if (Fortran_COMPILER_NAME MATCHES "ifort.*")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost -mkl -diag-disable 13003")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
# elseif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
# else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#   message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#   message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#   message ("No optimized Fortran compiler flags are known, we just try -O2...")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
# endif (Fortran_COMPILER_NAME MATCHES "ifort.*")

# build executables

# set (NMPROGRAMS "diagstd" "diagrtb" "proj_modes_bin" "pdbmat")
# set (EXECUTABLES "NORMA.exe" ${NMPROGRAMS})
# set (SCRIPTS "gen_pert.sh" "pert_multi_mode.sh")
#
# add_executable ("NORMA.exe" "NORMA.f")
# foreach (p ${NMPROGRAMS})
#   add_executable (${p} "${p}.f")
# endforeach (p)
#
# # install executables and scripts
# install (TARGETS ${EXECUTABLES}
#          RUNTIME DESTINATION "bin")
# install (PROGRAMS ${SCRIPTS}
#          DESTINATION "bin")
#-------------------------------------------------------------------------------
