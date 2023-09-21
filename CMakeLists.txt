#Projectname
project(polyde-fem Fortran)

#Minimum required Cmake Version
cmake_minimum_required(VERSION 3.24)

#Enables Fortran support in Cmake
enable_language(Fortran)
set(CMAKE_C_COMPILER /opt/intel/oneapi/compiler/2023.1.0/linux/bin/icx)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})

#  W I N D O W S
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("Fortran compiler ID: " ${CMAKE_Fortran_COMPILER_ID})
# Intel Composer
    IF(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel"))
# Sets the flags for all Build types Fortran Compiler
        SET(CMAKE_Fortran_FLAGS "/nologo /MP /libs:static /Qmkl:parallel")
# Sets the linker flags (Stack: 3072MB)
        SET(CMAKE_EXE_LINKER_FLAGS "/NOLOGO /INCREMENTAL:NO /STACK:3221225472")
# Sets the DEBUG flags for the Fortran Compiler
        SET(CMAKE_Fortran_FLAGS_DEBUG "/Od /fpe:0 /debug:full /traceback /dbglibs")
        SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "/DEBUG")
# Sets the MINSIZEREL flags for the Fortran Compiler (Threaded Debug)
        SET(CMAKE_Fortran_FLAGS_MINSIZEREL "/Od /fpe:0 /debug:full /traceback /Qopenmp /dbglibs /threads")
        SET(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "/DEBUG")
# Sets the RELWITHDEBINFO flags for the Fortran Compiler
        SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "/O2 /Ob2 /debug:full /Qopenmp /QxHost /Qip /assume:buffered_io /threads")
        SET(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "/DEBUG")
# Sets the RELEASE flags for the Fortran Compiler
        SET(CMAKE_Fortran_FLAGS_RELEASE "/O2 /Ob2 /Qopenmp /QxHost /Qipo /align:all /assume:buffered_io /threads")
# PGI Compiler
    ELSEIF(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI"))
# Flang Compiler
    ELSEIF(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang"))
# GFortran (GNU GCC) Compiler
    ELSEIF(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU"))
    INCLUDE_DIRECTORIES(
          ${CMAKE_BINARY_DIR}/Source/Types_Interfaces
          ${CMAKE_BINARY_DIR}/Packages/FunctionParser
          ${CMAKE_BINARY_DIR}/Packages/Json
          ${CMAKE_BINARY_DIR}/Source/Femlib
          ${CMAKE_BINARY_DIR}/Source/Femlib2D
          ${CMAKE_BINARY_DIR}/Source/Femlib3D
          ${CMAKE_BINARY_DIR}/Source/Post3D
          ${CMAKE_BINARY_DIR}/Source/System
          )
# Sets the flags for all Build types Fortran Compiler
        SET(CMAKE_Fortran_FLAGS "-static-libgfortran -ffixed-line-length-none -ffree-line-length-none -fdec-math")
# Sets the linker flags (Stack: 3072MB)
# -frecursive    : Allow indirect recursion by forcing all local arrays to be allocated on the stack.
#                  This flag cannot be used together with -fmax-stack-var-size= or -fno-automatic.
# -fstack-arrays : Adding this option will make the Fortran compiler put all arrays of unknown size and array temporaries onto stack memory.
        SET(CMAKE_EXE_LINKER_FLAGS "-fmax-stack-var-size=3221225472")
# Sets the DEBUG flags for the Fortran Compiler
        SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -warn all")
        #SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "/DEBUG")
# Sets the MINSIZEREL flags for the Fortran Compiler (Threaded Debug)
# To really represent a Min-Size-Release, this configuration should have optimisation level /-Os/, which
# optimizes similar to /-O2/ leaving out those optimizations that increase the executable size
#        SET(CMAKE_Fortran_FLAGS_MINSIZEREL "-g -Os -fcheck=all -fopenmp -fbacktrace /Qopenmp /dbglibs /threads")
# -ffpe-summary=all : List floating-point exceptions
#                     exceptions: ‘invalid’, ‘zero’, ‘overflow’, ‘underflow’, ‘inexact’ and ‘denormal’
#                     By default, a summary for all exceptions but ‘inexact’ is shown
        SET(CMAKE_Fortran_FLAGS_MINSIZEREL "-g -O0 -fcheck=all -fopenmp -fbacktrace")
        #SET(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "/DEBUG")
# Sets the RELWITHDEBINFO flags for the Fortran Compiler
# -fdump-fortran-optimized : Output the parse tree after front-end optimization.
        SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -fopenmp -fbacktrace")
        #SET(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "/DEBUG")
# Sets the RELEASE flags for the Fortran Compiler
# -fcoarray=lib     : Library-based coarray parallelization; a suitable GNU Fortran coarray library needs to be linked.
# -frepack-arrays   : In some circumstances GNU Fortran may pass assumed shape array sections via a descriptor describing a noncontiguous area of memory.
#                     This option adds code to the function prologue to repack the data into a contiguous block at runtime.
#                     This should result in faster accesses to the array. However it can introduce significant overhead to the function call, especially when the passed data is noncontiguous.
# -msse2, -mavx512f : Usage of specific ISA vectorization extensions (-march=native use the ones that are native to the compiling machine)
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native")
    ENDIF(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel"))
    message("${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} ${Fortran_COMPILER_NAME}")

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


#  L I N U X
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# For oneAPI MKL
  include_directories(/opt/intel/oneapi/mkl/latest/include)
# OneMKL setup with Intel OneAPI
# using config mode
# Looks for MKLConfig.cmake in the path
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
message("MKL variables from the file: " ${MKL_RPATH})

# PETSc setup with FindPkgConfig
# As in the sample in <petsc root>/petsc/share/petsc/CMakeLists.txt
#  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# set root of location to find PETSc's pkg-config
#set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
#set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)
#find_package(PkgConfig REQUIRED)
# part of PkgConfig
#pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
message("petsc variables: " ${PETSC})
# End PETSc setup

  include_directories(
    ${CMAKE_BINARY_DIR}/Source/Types_Interfaces
    ${CMAKE_BINARY_DIR}/Packages/FunctionParser
    ${CMAKE_BINARY_DIR}/Packages/Json
    ${CMAKE_BINARY_DIR}/Source/Femlib
    ${CMAKE_BINARY_DIR}/Source/Femlib2D
    ${CMAKE_BINARY_DIR}/Source/Femlib3D
    ${CMAKE_BINARY_DIR}/Source/Post3D
    ${CMAKE_BINARY_DIR}/Source/System
  )
# FFLAGS depend on the compiler
    message("Determining Fortran compiler:")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("Fortran compiler ID: " ${CMAKE_Fortran_COMPILER_ID})
    message ("Fortran compiler version: " ${CMAKE_Fortran_COMPILER_VERSION})

    if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
# gfortran
        set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
        set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
    elseif (Fortran_COMPILER_NAME MATCHES "ifx")
# ifort
# option -ip not applicable to ifx
# Sets the DEBUG flags for the Fortran Compiler
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fpe0 -debug full -traceback -qmkl -extend-source -check pointer -check uninit")
# Sets the MINSIZEREL flags for the Fortran Compiler (Threaded Debug)
        SET(CMAKE_Fortran_FLAGS_MINSIZEREL "-O0 -fpe0 -debug full -qmkl -qopenmp -extend-source -check pointer -check uninit -traceback")
# Sets the RELWITHDEBINFO flags for the Fortran Compiler
        SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -debug full -qmkl -qopenmp -xHost -inline-level=2 -extend-source -check pointer -check uninit")
# Sets the RELEASE flags for the Fortran Compiler
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -qmkl -qopenmp -xHost -inline-level=2 -extend-source -real-size 32")
# PolyDE-MTSO-specific Release compiler flags to meet the target systems needs (Intel(R) Xeon(R) CPU E5-2630 v4)
#        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -qmkl -qopenmp -axcore-avx2 -inline-level=2 -extend-source -real-size 32")
    elseif (Fortran_COMPILER_NAME MATCHES "g77")
# g77
        set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
        set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
    elseif (Fortran_COMPILER_NAME MATCHES "f95")
# f95
# Sets the flags for all Build types Fortran Compiler
# For explanations on the settings, please refer to the windows setup as provided above!
        SET(CMAKE_Fortran_FLAGS "-static-libgfortran -ffixed-line-length-none -ffree-line-length-none -fdec-math")
        SET(CMAKE_EXE_LINKER_FLAGS "-fmax-stack-var-size=3221225472")
        SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -Wall -Wextra -Wconversion-extra -Warray-temporaries -Wuse-without-only")
        SET(CMAKE_Fortran_FLAGS_MINSIZEREL "-g -O0 -fcheck=all -fopenmp -fbacktrace")
        SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -fopenmp -fbacktrace")
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native")
    else (Fortran_COMPILER_NAME MATCHES "gfortran")
# gfortran
# Sets the flags for all Build types Fortran Compiler
# For explanations on the settings, please refer to the windows setup as provided above!
        SET(CMAKE_Fortran_FLAGS "-static-libgfortran -ffixed-line-length-none -ffree-line-length-none -fdec-math")
        SET(CMAKE_EXE_LINKER_FLAGS "-fmax-stack-var-size=3221225472")
        SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -Wall -Wextra -Wconversion-extra -Warray-temporaries -Wuse-without-only")
        SET(CMAKE_Fortran_FLAGS_MINSIZEREL "-g -O0 -fcheck=all -fopenmp -fbacktrace")
        SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -fopenmp -fbacktrace")
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native")
    endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# Sets the build type to Release/Debug
# If uncommented, this is set by the cmake builder (e.g. VSCode cmake plugin)
set(CMAKE_BUILD_TYPE "Release") # "Debug" or "Release". See other options
message ("Active configuration: " ${CMAKE_BUILD_TYPE})

#Subdirectories in which other CMakeLists.txt or code files are located
add_subdirectory(Packages)
add_subdirectory(Source)
