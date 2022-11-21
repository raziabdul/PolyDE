#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)
#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib3D)

#Add an executable to the project using the source files from SRC 
add_executable(post3d ${SRC})

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(post3d Femlib Femlib3D System)

#Link a target to given libraries 
#target_link_libraries(Post3D Femlib Femlib3D Pgplot System)
target_link_libraries(post3d Femlib Femlib3D pgplot System X11 png)

#Additional compile flags for specifying the runtime library 
#set_target_properties(Post3D PROPERTIES COMPILE_FLAGS "/warn:all /libs:static /threads /dbglibs")
set_target_properties(post3d PROPERTIES COMPILE_FLAGS "-warn all -threads")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS post3d
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)

