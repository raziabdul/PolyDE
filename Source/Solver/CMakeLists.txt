# Disable pardiso_unsym.f90 for now
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
	preassemb.f90
	residual.f90
	solve2.f90
	solver.f90
	transientsolver.f90
	umfsolver.f90
)
#Additional files for the current project are stored in variable SRC2
SET(SRC2 
	${CMAKE_SOURCE_DIR}/Source/Post/zanfpp.f90
	${CMAKE_SOURCE_DIR}/Source/Post/circle.f90
	${CMAKE_SOURCE_DIR}/Source/Post/farbe.f90
	${CMAKE_SOURCE_DIR}/Source/Post/farbpalette.f90
	${CMAKE_SOURCE_DIR}/Source/Post/hsb_to_rgb.f90
	${CMAKE_SOURCE_DIR}/Source/Post/struc.f90
	${CMAKE_SOURCE_DIR}/Source/Post/triplt.f90
)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)

#Add an executable to the project using the source files from SRC and SRC2
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_executable(Solver WIN32 ${SRC} ${SRC2})
else()
  add_executable(solver ${SRC} ${SRC2})
endif()

# Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(solver Femlib Fluidlib System)

# Link a target to given libraries 
# On Linux, Blas must appear after Umfpack in the list
target_link_libraries(solver Femlib pgplot Fluidlib Umfpack System 
	X11 png PkgConfig::PETSC)

#Additional compile flags for using the "Intel Math Kernel Library" 
#and for specifying the runtime library 
#set_target_properties(Solver PROPERTIES COMPILE_FLAGS "/warn:all /Qmkl:parallel /libs:qwin /dbglibs")
set_target_properties(solver PROPERTIES COMPILE_FLAGS "-warn all")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS solver
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)