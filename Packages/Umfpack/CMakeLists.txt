#Generate a list of files with filetype .f or .f90 which are in 
#the current directory and store it in the variable SRC
file(GLOB SRC *.f *.f90 *.for)

#Add a library to the project using the source files from SRC 
add_library(Umfpack ${SRC})

#Additional compile flags for specifying the runtime library 
#set_target_properties(Umfpack PROPERTIES COMPILE_FLAGS "/libs:static /threads /dbglibs")
set_target_properties(Umfpack PROPERTIES COMPILE_FLAGS "-threads")

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS Umfpack
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)

