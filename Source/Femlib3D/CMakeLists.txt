#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90)

#Additional files for the current project are stored in variable SRC2
SET(SRC2 
	${CMAKE_SOURCE_DIR}/Source/Userinterface3D/userbc3D.f90
	${CMAKE_SOURCE_DIR}/Source/Userinterface3D/usermaterials3D.f90
)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/Polyde/build/Source/Femlib)
include_directories(${CMAKE_SOURCE_DIR}/Polyde/build/Packages/Json)

#Add a library to the project using the source files from SRC and SRC2
add_library(Femlib3D ${SRC} ${SRC2})

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(Femlib3D Femlib Json)

#Link a target to given libraries 
target_link_libraries(Femlib3D Femlib Json)

#Additional compile flags for specifying the runtime library 
#set_target_properties(Femlib3D PROPERTIES COMPILE_FLAGS "/warn:all /libdir:noauto /libs:static /threads /dbglibs")
set_target_properties(Femlib3D PROPERTIES COMPILE_FLAGS "-warn all -threads")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS Femlib3D
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
