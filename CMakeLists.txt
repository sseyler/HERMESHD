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
# set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} cmake")
# include(ConfigurePythonWrapping)
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
#-------------------------------------------------------------------------------


################################################################################
# default installation
#********************************************
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
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
set(default_sources ${sources})

# NOTE: It will be useful to have separate subdirectories, probably
list(REMOVE_ITEM default_sources ${SOURCE_DIR}/hermeshd.f90)
list(REMOVE_ITEM default_sources ${SOURCE_DIR}/f90wrap_hermeshd.f90)
list(REMOVE_ITEM default_sources ${SOURCE_DIR}/spatial.f90)

# set (EXECUTABLE "hermes")
# add_executable(${EXECUTABLE} ${default_sources})

# install (TARGETS ${EXECUTABLES} RUNTIME DESTINATION "bin") # already handled by EXECUTABLE_OUTPUT_PATH
#-------------------------------------------------------------------------------


find_package(PythonInterp REQUIRED)
find_program(F90WRAP_EXEC NAMES "f90wrap" REQUIRED)
find_program(F2PY_F90WRAP_EXEC NAMES "f2py-f90wrap" REQUIRED)

set(module_sources ${sources})
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/main.f90)
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/spatial.f90)
list(REMOVE_ITEM module_sources ${SOURCE_DIR}/f90wrap_hermeshd.f90)


set (PYFILES "hermeshd.py")


# set(F2PY_ADD_LIBS "mkl_intel_lp64 mkl_sequential mkl_core mkl_vml_avx")
# message("F2PY_ADD_LIBS:  ${F2PY_ADD_LIBS}")

set(MKLROOT /nfs/packages/opt/Linux_x86_64/intel/13.0/mkl)
set(MKL_TARGET_ARCH intel64)
set(MKLPATH ${MKLROOT}/lib/${MKL_TARGET_ARCH})


function(F2PY_ADD)

    include(CMakeParseArguments)

    cmake_parse_arguments(F2PY_ADD "DONT_COMPILE_SRC;SIGNATURE" "MODULE;DESTINATION" "SRC;LIBS;EXTERNAL_LIBS" ${ARGN})

    # check the arguements:
    # - check that MODULE is defined
    # - check that SRC is defined
    # - check that if DONT_COMPILE_SRC and SIGNATURE is defined that SRC
    if (NOT F2PY_ADD_MODULE)
        message(FATAL_ERROR "[F2PY] Undefined module in f2py_add. Add the module name to the arguements f2py_add(MODULE foo ..)")
    endif()
    if (NOT F2PY_ADD_SRC)
        message(FATAL_ERROR "[F2PY] Undefined source in f2py_add. Add the source to the arguments f2py_add(... SRC foo.f90..)")
    endif()
    if (F2PY_ADD_DONT_COMPILE_SRC AND NOT F2PY_ADD_LIBS)
        message(WARNING "[F2PY] f2py_add is called with DONT_COMPILE_SRC but with out any LIBS")
    endif()

    # give information
    message(STATUS "Building f2py module ${F2PY_ADD_MODULE} with ${F2PY_ADD_SRC}")
    if (F2PY_ADD_EXTERNAL_LIBS)
        message(STATUS "F2PY: Building f2py module ${F2PY_ADD_MODULE} with external lib ${F2PY_ADD_EXTERNAL_LIBS}")
    endif()
    if (F2PY_ADD_LIBS)
        message(STATUS "F2PY: Building f2py module ${F2PY_ADD_MODULE} with libraries ${F2PY_ADD_LIBS}")
    endif()
    if (F2PY_ADD_SIGNATURE)
        message(STATUS "F2PY: Building f2py module ${F2PY_ADD_MODULE} with custom signature")
    endif()


    # Make the source filenames absolute. The result is in _abs_srcs
    set(_abs_srcs)
    foreach(_src ${module_sources})
        get_filename_component(_abs_src ${_src} ABSOLUTE)
        list(APPEND _abs_srcs ${_abs_src})
    endforeach()

    # Make a line like -llib1 -l1ib2 -llib3 etc. The result is in _flibs
    if (F2PY_ADD_LIBS)
        set(_flibs)
        foreach(_lib ${F2PY_ADD_LIBS})
            list(APPEND _flibs -l${_lib})
            get_target_property(_liblocation ${_lib} LOCATION)
            get_filename_component(_libpath ${_liblocation} PATH)
            find_library(_libpath ${_lib})
            message("Found the library path ${_liblocation}")
            message("Found the path ${_libpath}")
        endforeach()
    endif()


    add_custom_command(
        OUTPUT  f90wrap_${F2PY_ADD_MODULE}.f90
        COMMAND ${F90WRAP_EXEC} -m ${F2PY_ADD_MODULE} ${SOURCE_DIR}/${F2PY_ADD_MODULE}.f90
                --only main setup step cleanup generate_output
        DEPENDS ${SOURCE_DIR}/${F2PY_ADD_MODULE}.f90
        COMMENT "${F90WRAP_EXEC}: wrapping source files and generating Python module..."
    )


    # add_custom_command(
    #     OUTPUT  gen
    #     COMMAND ${CMAKE_Fortran_COMPILER} -c -fPIC ${module_sources} > gen
    #     DEPENDS ${SOURCE_DIR}/${F2PY_ADD_MODULE}.f90
    #     COMMENT "Compiling Fortran modules and generating object files..."
    # )

    # file(GLOB mod_files ${CMAKE_CURRENT_SOURCE_DIR}/*.mod)
    set(_depends f90wrap_${F2PY_ADD_MODULE}.f90)

    add_custom_command(
        OUTPUT  _${F2PY_ADD_MODULE}.so
        COMMAND ${F2PY_F90WRAP_EXEC} -m _${F2PY_ADD_MODULE}
                --f90exec=${CMAKE_Fortran_COMPILER}
                --opt="-O2"	--f90flags=${CMAKE_Fortran_FLAGS}
                -c f90wrap_${F2PY_ADD_MODULE}.f90 -L${CMAKE_CURRENT_BINARY_DIR} -l${F2PY_ADD_MODULE}
                -L${MKLPATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx
        DEPENDS ${_depends}
        COMMENT "${F2PY_F90WRAP_EXEC}: compiling extension modules..."
    )


    add_custom_target(_${F2PY_ADD_MODULE} ALL DEPENDS _${F2PY_ADD_MODULE}.so)

endfunction(F2PY_ADD)



set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${flags} -fPIC")
add_library(hermeshd ${module_sources})

F2PY_ADD(
    MODULE hermeshd
    SRC ${module_sources}
    EXTERNAL_LIBS "-L${MKLPATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx"
)


add_dependencies(_hermeshd hermeshd)


# install (TARGETS _${F2PY_ADD_MODULE}.so RUNTIME DESTINATION ${SCRATCH})
# install (TARGETS ${PYFILES} RUNTIME DESTINATION ${SCRATCH})
