################################################################################
# CMake include file for wrapping HERMESHD with Python
################################################################################

################################################################################
# Set Fortran compiler for f2py
#********************************************
# macro (set_fcompiler)
#     if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#         set(F2PY_FCOMPILER "intelem")
#     endif()
# endmacro (set_fcompiler)
#-------------------------------------------------------------------------------


################################################################################
# Wrap Fortran in Python
#********************************************
# Usage
#--------------------------------------------
# PYWRAP (
#     MOD <module_name>
#     SRC <list of sources>
#     SUBPROGRAMS <list of subprograms to wrap>
#     LIBS <local libraries to link>
#     EXT_LIBS <external libraries to link>
# )
#
# SRC takes a list of source files: "srcfile1.f90 srcfile2.f90 ..."
# SUBPROGRAMS takes a list of subprogram names to be wrapped w/ f90wrap: "subprogram1 subprogram2 ..."
# LIBS takes a list of library names (w/o "lib" preprended): "nameoflib1 nameoflib2 ..."
# EXT_LIBS should take the explicit form: "-L<path/to/lib> -lextlib1 -lextlib2 ..."
#-------------------------------------------------------------------------------
macro(PYWRAP)

    include(CMakeParseArguments)
    cmake_parse_arguments(PYWRAP "NO_COMPILE" "MOD" "SRC;SUBPROGRAMS;LIBS;EXT_LIBS" ${ARGN})

    find_package(PythonInterp REQUIRED)  # Find the Python interpreter
    find_program(F90WRAP_EXEC NAMES "f90wrap" REQUIRED)  # Get the f90wrap executable
    find_program(F2PY_F90WRAP_EXEC NAMES "f2py-f90wrap" REQUIRED)  # Get the f2py-f90wrap executable

    #***********************************************************
    # Check for inputs
    #********************************************
    if (NOT PYWRAP_MOD)
        message(FATAL_ERROR "MOD is undefined in PYWRAP: add the module name to the arguments using PYWRAP(MOD <name> ...)")
    endif()
    if (NOT PYWRAP_SRC)
        message(FATAL_ERROR "SRC is undefined in PYWRAP: add source(s) to the arguments using PYWRAP(... SRC <source> ...)")
    endif()
    if (PYWRAP_NO_COMPILE AND NOT PYWRAP_LIBS)
        message(WARNING "NO_COMPILE specified in PYWRAP: LIBS is undefined")
    endif()
    #-----------------------------------------------------------


    #***********************************************************
    # Print out information
    #********************************************
    message(STATUS "Building module ${PYWRAP_MOD} from the following sources: ${PYWRAP_SRC}")
    if (PYWRAP_LIBS)
        message(STATUS " >>> ${PYWRAP_MOD} will be built from the following libraries: ${PYWRAP_LIBS}")
    endif()
    if (PYWRAP_EXT_LIBS)
        message(STATUS " >>> ${PYWRAP_MOD} will use the following external libraries: ${PYWRAP_EXT_LIBS}")
    endif()
    #-----------------------------------------------------------


    #***********************************************************
    # Generate wrapper (as input to f2py/f2py-f90wrap)
    #********************************************
    # NOTE: also generates ${PYWRAP_MOD}.py
    #--------------------------------------------
    add_custom_command(
        OUTPUT  f90wrap_${PYWRAP_MOD}.f90
        COMMAND ${F90WRAP_EXEC} -m ${PYWRAP_MOD} ${SOURCE_DIR}/${PYWRAP_MOD}.f90
                --only ${PYWRAP_SUBPROGRAMS}
        DEPENDS ${SOURCE_DIR}/${PYWRAP_MOD}.f90
        COMMENT "${F90WRAP_EXEC}: wrapping source files and generating Python module..."
    )
    #-----------------------------------------------------------


    #***********************************************************
    # Compile extension module
    #********************************************
    if (NOT PYWRAP_LIBS)
        set(libs_list "-l${PYWRAP_MOD}")
    # else (PYWRAP_LIBS)
    #     set(libs_list)
    #     foreach(lib ${PYWRAP_LIBS})
    #         list(APPEND libs_list -l${lib})
    #         get_target_property(lib_loc ${lib} LOCATION)
    #         get_filename_component(lib_path ${lib_loc} PATH)
    #         find_library(lib_path ${lib})
    #         message("Found the library path ${lib_loc}")
    #         message("Found the path ${lib_path}")
    #     endforeach()
    endif()

    # file(GLOB mod_files ${CMAKE_CURRENT_SOURCE_DIR}/*.mod)
    if (PYWRAP_NO_COMPILE)
        set(_depends f90wrap_${PYWRAP_MOD}.f90 )
    else()
        set(_depends f90wrap_${PYWRAP_MOD}.f90)
    endif()

    add_custom_command(
        OUTPUT  _${PYWRAP_MOD}.so
        COMMAND ${F2PY_F90WRAP_EXEC} -m _${PYWRAP_MOD}
                --f90exec=${CMAKE_Fortran_COMPILER}
                --opt="-O2"	--f90flags=${CMAKE_Fortran_FLAGS}
                -c f90wrap_${PYWRAP_MOD}.f90 -L${CMAKE_CURRENT_BINARY_DIR} ${libs_list}
                ${PYWRAP_EXT_LIBS}
        DEPENDS ${_depends}
        COMMENT "${F2PY_F90WRAP_EXEC}: compiling extension modules..."
    )
    #-----------------------------------------------------------

    add_custom_target(_${PYWRAP_MOD} ALL DEPENDS _${PYWRAP_MOD}.so)

endmacro(PYWRAP)
#-------------------------------------------------------------------------------
