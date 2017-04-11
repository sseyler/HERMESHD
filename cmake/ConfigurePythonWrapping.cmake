################################################################################
# CMake include file for wrapping HERMESHD with Python
################################################################################


find_package(PythonInterp REQUIRED)


## Get the f90wrap and f2py-f90wrap executables

find_program(F90WRAP_EXEC NAMES "f90wrap" REQUIRED)
find_program(F2PY_F90WRAP_EXEC NAMES "f2py-f90wrap" REQUIRED)


################################################################################
# Set Fortran compiler for f2py
#********************************************
macro (set_fcompiler)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(F2PY_FCOMPILER "intelem")
    endif()
endmacro (set_fcompiler)
#-------------------------------------------------------------------------------


################################################################################
# Generate wrapper (as input to f2py/f2py-f90wrap)
#********************************************
macro (generate_wrapper mod_name)
    add_custom_command(
        OUTPUT "${mod_name}"
        COMMAND ${F90WRAP_EXEC} -m ${mod_name}
                ${sources} --only ${pymodsubprograms}
        DEPENDS ${sources}
        COMMENT "${F90WRAP_EXEC}: wrapping source files and generating Python module"
    )
endmacro (generate_wrapper)
#-------------------------------------------------------------------------------


################################################################################
# Compile modules and generate object files
#********************************************
macro (compile_modules _sources)

    file(GLOB sources src/*.f90)
    set(f2py_sources ${sources})
    list(REMOVE_ITEM f2py_sources ${CMAKE_CURRENT_SOURCE_DIR}/main.f90)

    add_custom_command(
        OUTPUT "${mod_name}"
        COMMAND ${CMAKE_Fortran_COMPILER} -c -fPIC
                 ${CMAKE_Fortran_FLAGS} ${_sources}
        DEPENDS ${_sources}
        COMMENT "${F90WRAP_EXEC}: wrapping source files and generating Python module"
    )
endmacro (generate_wrapper)
#-------------------------------------------------------------------------------


################################################################################
# Compile extension module
#********************************************
macro (compile_extension _sources)

    set(MKLROOT /nfs/packages/opt/Linux_x86_64/intel/13.0/mkl)
    set(MKL_TARGET_ARCH intel64)
    set(MKLPATH ${MKLROOT}/lib/${MKL_TARGET_ARCH})

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} ${flags} -fc=ifort")

    file(GLOB f90wrap_sources src/f90wrap_*.f90)
    file(GLOB objects ${CMAKE_Fortran_MODULE_DIRECTORY}/*.o)

    add_custom_command(
        OUTPUT "${SOME_KIND_OF_NAME}"
        COMMAND ${F2PY_F90WRAP_EXEC} -m _${mod_name}
                --build-dir "${CMAKE_CURRENT_SOURCE_DIR}"
                --f90exec=${CMAKE_Fortran_COMPILER}
                --opt="-O2"	--f90flags=${CMAKE_Fortran_FLAGS}
                -c ${objects} ${f90wrap_sources}
                -L${MKLPATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_vml_avx
        DEPENDS ${sources}
        COMMENT "${F2PY_F90WRAP_EXEC}: compiling extension modules"
    )
endmacro (compile_extension)
#-------------------------------------------------------------------------------


file(GLOB sources src/*.f90)
set(module_sources ${sources})
list(REMOVE_ITEM module_sources ${CMAKE_CURRENT_SOURCE_DIR}/main.f90)
