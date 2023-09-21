#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90)

#Additional files for the current project are stored in variable SRC2
SET(SRC2
	${CMAKE_SOURCE_DIR}/Source/Solver/getbcval2D.f90
	${CMAKE_SOURCE_DIR}/Source/Solver/residual.f90
)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)

#Add an executable to the project using the source files from SRC and SRC2
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_executable(Post WIN32 ${SRC} ${SRC2})
else()
  add_executable(post ${SRC} ${SRC2})
endif()

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
# Remove Pgplot. Not needed on Linux since pre-installed on /usr/lib64 
add_dependencies(post Femlib System Userinterface)

#Link a target to given libraries 
#target_link_libraries(Post Femlib pgplot System Userinterface)
target_link_libraries(post Femlib pgplot System Userinterface X11 png)

#Additional compile flags for specifying the runtime library 
#set_target_properties(Post PROPERTIES COMPILE_FLAGS "/warn:all /libs:qwin")
set_target_properties(post PROPERTIES COMPILE_FLAGS "-warn all")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS post
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
