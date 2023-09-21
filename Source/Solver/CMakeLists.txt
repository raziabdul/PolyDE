# Include manually list of files for solver
SET(SRC
	adaptation.f90
	assemblyandsolve.f90
	assemblyandsolvetrans.f90
	assembly.f90
	assemblytrans.f90
	coloring.f90
	csrmmout.f90
	elementmatrix.f90
	elementmatrixtrans.f90
	getbcval2D.f90
	hadapt.f90
	hpadapt.f90
	hpffein.f90
	itrhstout.f90
	matrhsout.f90
	mmio.f
	nonlinearsolver.f90
	padapt.f90
	pardiso_unsym.f90
	preassemb.f90
	residual.f90
	solve2.f90
	solver.f90
	transientsolver.f90
) 

# Additional files for "solver" are stored in variable SRC2
SET(SRC2 
	${CMAKE_SOURCE_DIR}/Source/Post/zanfpp.f90
	${CMAKE_SOURCE_DIR}/Source/Post/circle.f90
	${CMAKE_SOURCE_DIR}/Source/Post/farbe.f90
	${CMAKE_SOURCE_DIR}/Source/Post/farbpalette.f90
	${CMAKE_SOURCE_DIR}/Source/Post/hsb_to_rgb.f90
	${CMAKE_SOURCE_DIR}/Source/Post/struc.f90
	${CMAKE_SOURCE_DIR}/Source/Post/triplt.f90
)

message("Check MKL_INCLUDE: " ${MKL_INCLUDE})

# Must have. Use add_library to build the MKL Pardiso interface as a module
add_library(PardisoMod ${MKL_INCLUDE}/mkl_pardiso.f90)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)

#Add an executable to the project using the source files from SRC and SRC2
add_executable(solver ${SRC} ${SRC2})

# Add a dependency between targets (target-name depend-target1 depend-target2 ...)
# Remove Pgplot, Qmrpack, UmfPack, Blas dependencies to make build faster
# add_dependencies(Solver Femlib Blas  Fluidlib  Umfpack Umfpackclib System)
add_dependencies(solver Femlib Fluidlib System)

# pontentially necessary
#libmkl_blas95_lp64.a
#libmkl_core.a
#libmkl_intel_lp64.a
#libmkl_intel_thread.a
#libmkl_lapack95_lp64.a
#libmkl_scalapack_lp64.a
#libmkl_sequential.a
#libmkl_sycl.a
#libmkl_tbb_thread.a

# Link a target to given libraries 
# Tested with: mkl_core mkl_intel_lp64 mkl_intel_thread mkl_lapack95_lp64 mkl_blas95_lp64 mkl_gf_lp64 mkl_sequential mkl_sycl mkl_tbb_thread

target_link_libraries(solver Femlib pgplot Fluidlib System 
	X11 png  mkl_core)
# to enable petsc add this above:	PkgConfig::PETSC
# Additional compile flags for using the "Intel Math Kernel Library" 
# and for specifying the runtime library 
# set_target_properties(Solver PROPERTIES COMPILE_FLAGS "/warn:all /Qmkl:parallel /libs:qwin /dbglibs")
set_target_properties(solver PROPERTIES COMPILE_FLAGS "-warn all")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS solver
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)

