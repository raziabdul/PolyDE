#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)
#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib3D)

#Add a library to the project using the source files from SRC 
add_library(Userinterface3D ${SRC})

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(Userinterface3D Femlib Femlib3D)

#Link a target to given libraries 
target_link_libraries(Userinterface3D Femlib Femlib3D)

#Additional compile flags for specifying the runtime library 
#set_target_properties(Userinterface3D PROPERTIES COMPILE_FLAGS "/warn:all /libdir:noauto /libs:static /threads /dbglibs")
set_target_properties(Userinterface3D PROPERTIES COMPILE_FLAGS "-warn all -threads")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS Userinterface3D
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
