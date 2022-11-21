#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90)

#Additional files for the current project are stored in variable SRC2
SET(SRC2 
	${CMAKE_SOURCE_DIR}/Source/Post/circle.f90
	${CMAKE_SOURCE_DIR}/Source/Post/farbe.f90
	${CMAKE_SOURCE_DIR}/Source/Post/farbpalette.f90
	${CMAKE_SOURCE_DIR}/Source/Post/hsb_to_rgb.f90
	${CMAKE_SOURCE_DIR}/Source/Post/struc.f90
	${CMAKE_SOURCE_DIR}/Source/Post/triplt.f90
	${CMAKE_SOURCE_DIR}/Source/Post/zanfpp.f90
)

#Add path to modules
include_directories(${CMAKE_SOURCE_DIR}/build/Source/Femlib)

#Add an executable to the project using the source files from SRC and SRC2
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_executable(Triang WIN32 ${SRC} ${SRC2})
else()
  add_executable(triang ${SRC} ${SRC2})
endif()

#Add a dependency between targets (target-name depend-target1 depend-target2 ...)
add_dependencies(triang Femlib System)

#Link a target to given libraries 
target_link_libraries(triang Femlib pgplot System X11 png)

#Additional compile flags for specifying the runtime library 
#set_target_properties(Triang PROPERTIES COMPILE_FLAGS "/warn:all /libs:qwin")
set_target_properties(triang PROPERTIES COMPILE_FLAGS "-warn all")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS triang
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)