#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
SET(SRC
	${CMAKE_SOURCE_DIR}/Source/Solver3D/assembly3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/coloring3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/CSR_residual.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/element_mark3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/elementmatrix3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/elementmatrix_scalar3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/estimate_error3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/getbcval.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/getepbc.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/getkpv.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/get_modres.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/get_pconfirm.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/h_adapt3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/hmesh_refine3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/hp_adapt3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/hp_alg_F8R.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/p_adapt3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/pmesh_refine3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/preassemb3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/prepare_adaptation3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/residual3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/residualsc3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/solve_adapt3D.F90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/solve_linear3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/solve_nonlinear3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/solver3D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver3D/vtk_varmatplot.f90
	)
#	${CMAKE_SOURCE_DIR}/Source/Solver3D/petscNonlinear.F90

#Additional files for the current project are stored in variable SRC2
SET(SRC2
	${CMAKE_SOURCE_DIR}/Source/Solver/itrhstout.f90
	${CMAKE_SOURCE_DIR}/Source/Solver/matrhsout.f90
	${CMAKE_SOURCE_DIR}/Source/Solver/solve2.f90
	)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)
#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib3D)
# For Linux, OneMKL path already included in the main dir.

#Add an executable to the project using the source files from SRC and SRC2
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_executable(Solver3D WIN32 ${SRC} ${SRC2})
else()
  add_executable(solver3d ${SRC} ${SRC2})
endif()

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
# Remove Pgplot. Not needed on Linux since pre-installed on /usr/lib64 
# Remove Qmrpack dependency to make build faster

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	add_dependencies(Solver3D Femlib Femlib3D Geompack System)
else()
	add_dependencies(solver3d Femlib Femlib3D Geompack System)
endif()


#Link a target to given libraries 
#target_link_libraries(Solver3D Femlib Femlib3D Pgplot Qmrpack Umfpack Umfpackclib Geompack System)
# For petsc add to target: PkgConfig::PETSC
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(solver3d Femlib Femlib3D pgplot Geompack System X11)
else()
	message("To add")
endif()


#Additional compile flags for using the "Intel Math Kernel Library" 
#and for specifying the runtime library
# warn options: all, external 
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	set_target_properties(solver3d PROPERTIES COMPILE_FLAGS "-warn all -fpe3 -threads")
else()
	set_target_properties(Solver3D PROPERTIES COMPILE_FLAGS "/warn:all /fpe:0 /Qmkl:parallel /libs:static /threads /dbglibs")
endif()


#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	install(TARGETS solver3d
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)
else()
	message("To add")
endif()