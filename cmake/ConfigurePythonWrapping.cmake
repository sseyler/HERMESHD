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
# SUBPROGRAMS takes a string of names separated by spaces (subprograms to be wrapped w/ f90wrap): "subprogram1 subprogram2 ..."
# LIBS takes a list of library names (w/o "lib" preprended): nameoflib1 nameoflib2 ...
# EXT_LIBS takes a list of libraries: extlib1 extlib2 ...
# EXT_LIBS_LOC takes the path to the external libraries: <path/to/lib>
#-------------------------------------------------------------------------------
macro(PYWRAP)

    include(CMakeParseArguments)
    cmake_parse_arguments(PYWRAP "NO_COMPILE" "MOD" "SRC;SUBPROGRAMS;LIBS;EXT_LIBS;EXT_LIBS_LOC" ${ARGN})

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
    # Parse libraries and their locations
    #********************************************
    # if (NOT PYWRAP_LIBS)
    #     set(libs_list "-l${PYWRAP_MOD}")
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
    # endif()

    if (PYWRAP_EXT_LIBS)
        if (NOT PYWRAP_EXT_LIBS_LOC)
            message(WARNING "No path to external libraries specified")
        endif()
        set (ext_lib_loc ${PYWRAP_EXT_LIBS_LOC})
        set (ext_lib_list)
        foreach (ext_lib ${PYWRAP_EXT_LIBS})
            list (APPEND ext_lib_list -l${ext_lib})
            message("Found the external library location: ${ext_lib_loc}")
            message("Found the external library path: ${ext_lib}")
        endforeach()
    endif()
    #-----------------------------------------------------------


    #***********************************************************
    # Print out information
    #********************************************
    set (BARE_SRC_NAMES)
    foreach (PYWRAP_SRC_NAME ${PYWRAP_SRC})
        get_filename_component(BARE_PYWRAP_SRC ${PYWRAP_SRC_NAME} NAME)
        list (APPEND BARE_SRC_NAMES ${BARE_PYWRAP_SRC})
    endforeach()
    message(STATUS "Building module ${PYWRAP_MOD} from the following sources: \n     ${BARE_SRC_NAMES}")

    if (PYWRAP_LIBS)
        message(STATUS "Building from the following libraries: ${PYWRAP_LIBS}")
    endif()
    if (PYWRAP_EXT_LIBS AND PYWRAP_EXT_LIBS_LOC)
        message(STATUS "External libraries will be linked using: \n     -L${ext_lib_loc} ${ext_lib_list}")
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
    # Two choices (see below).
    # NOTE: The first (1) is preferable because it doesn't create an extra directory
    #       for the build. For (2), using add_custom_target to create "_${PYWRAP_MOD}"
    #       ends up generating a separate directory for it. However, using
    #       add_custom_command with "TARGET" allows the linking done by f2py-f90wrap
    #       to depend directly on the building of library ${PYWRAP_MOD} (so it
    #       ends up using the same directory ${PYWRAP_MOD}).

    # (1) EITHER THIS:
    #*************
    # add_custom_target(f90wrap_${PYWRAP_MOD} ALL DEPENDS f90wrap_${PYWRAP_MOD}.f90)
    # add_custom_command(
    #     TARGET  ${PYWRAP_LIBS} POST_BUILD
    #     COMMAND ${F2PY_F90WRAP_EXEC} -m _${PYWRAP_MOD}
    #             --f90exec=${CMAKE_Fortran_COMPILER}
    #             --opt="-O2"	--f90flags=${CMAKE_Fortran_FLAGS}
    #             -c f90wrap_${PYWRAP_MOD}.f90 -I${CMAKE_Fortran_MODULE_DIRECTORY}
    #             -L${CMAKE_CURRENT_BINARY_DIR} -l${PYWRAP_LIBS}  #${libs_list}
    #             -L${ext_lib_loc} ${ext_lib_list} #-L${MKLPATH} ${PYWRAP_EXT_LIBS}
    #     DEPENDS f90wrap_${PYWRAP_MOD}
    #     COMMENT "${F2PY_F90WRAP_EXEC}: compiling extension modules..."
    # )

    # (2) OR THIS:
    #*************
    if (PYWRAP_NO_COMPILE)
        set(_depends f90wrap_${PYWRAP_MOD}.f90 ${PYWRAP_MOD})
    else()
        set(_depends f90wrap_${PYWRAP_MOD}.f90 ${PYWRAP_MOD})
    endif()

    add_custom_command(
        OUTPUT  _${PYWRAP_MOD}.so
        COMMAND ${F2PY_F90WRAP_EXEC} --f90exec=${CMAKE_Fortran_COMPILER}
                --opt="-O2"	--f90flags=${CMAKE_Fortran_FLAGS}
                -c -m _${PYWRAP_MOD} f90wrap_${PYWRAP_MOD}.f90
                -I${CMAKE_Fortran_MODULE_DIRECTORY}
                -L${CMAKE_CURRENT_BINARY_DIR} -l${PYWRAP_LIBS}  #${libs_list}
                -L${ext_lib_loc} ${ext_lib_list} #-L${MKLPATH} ${PYWRAP_EXT_LIBS}
        DEPENDS ${_depends}
        COMMENT "${F2PY_F90WRAP_EXEC}: compiling extension modules..."
    )
    add_custom_target(_${PYWRAP_MOD} ALL DEPENDS _${PYWRAP_MOD}.so)
    add_dependencies(_${PYWRAP_MOD} ${PYWRAP_MOD})
    # Must do this if ${PYWRAP_MOD} isn't in ${_depends}:
    #   add_dependencies(_${PYWRAP_MOD} ${PYWRAP_MOD})
    # Can also make the custom target depend on the ${PYWRAP_MOD} library:
    #   add_custom_target(_${PYWRAP_MOD} ALL DEPENDS _${PYWRAP_MOD}.so ${PYWRAP_MOD})
    #-----------------------------------------------------------

    # install (FILES _${PYWRAP_MOD}.so DESTINATION ${SCRATCH}/lib)
    # install (TARGET ${PYWRAP_MOD} LIBRARY DESTINATION ${SCRATCH}/lib)
    # install (FILES ${PYWRAP_MOD}.py DESTINATION ${SCRATCH})

endmacro(PYWRAP)
#-------------------------------------------------------------------------------
