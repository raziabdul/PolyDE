#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)

#Add a library to the project using the source files from SRC 
add_library(Fluidlib ${SRC})

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(Fluidlib Femlib Userinterface)

#Link a target to given libraries 
target_link_libraries(Fluidlib Femlib Userinterface)

# Additional compile flags for specifying the runtime library 
set_target_properties(Fluidlib PROPERTIES COMPILE_FLAGS "-warn all -threads")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS Fluidlib
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)

