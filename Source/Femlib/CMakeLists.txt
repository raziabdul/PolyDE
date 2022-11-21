

# Generate a list of files with filetype .f or .f90 which are in 
# the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90 *.F90)

#Additional files for the current project are stored in variable SRC2
SET(SRC2 
	${CMAKE_SOURCE_DIR}/Source/Femlib3D/feminterface3d.f90
	${CMAKE_SOURCE_DIR}/Source/Post/fieldquantity.f90
	${CMAKE_SOURCE_DIR}/Source/Userinterface/usermaterials.f90
	${CMAKE_SOURCE_DIR}/Source/Userinterface3D/usermaterials3D.f90
    ${CMAKE_SOURCE_DIR}/Source/System/Linux/petsc/petscsolver.F90
)

# Subdirectories that contain PolyDE modules. 
# Petsc modules are taken care by cmake FindPkgConfig
include_directories(${CMAKE_SOURCE_DIR}/build/Packages/Json)
include_directories(${CMAKE_SOURCE_DIR}/build/Source/System)


# Add a library to the project using the source files from SRC and SRC2 
add_library(Femlib ${SRC} ${SRC2})

# Additional compile flags for specifying the runtime library 
# set_target_properties(Femlib PROPERTIES COMPILE_FLAGS "/warn:all /libdir:noauto /libs:static /threads /dbglibs")
set_target_properties(Femlib PROPERTIES COMPILE_FLAGS "-warn all -threads")

# Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(Femlib Json)

# Link a target to given libraries
message(${COMPILE_FLAGS} "  FEMLIB. Check PETsc")
target_link_libraries(Femlib Json PkgConfig::PETSC)

# Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS Femlib
	RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)