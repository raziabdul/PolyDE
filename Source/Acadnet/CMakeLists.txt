# Generate a list of files with filetype .f or .f90 which are in 
# the current directory and store it in the variable SRC

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  file(GLOB SRC *.f *.f90)
else()
  set(SRC
    acadnet.f90
    arc.f90
    blocks.f90
    circ.f90
    classes.f90
    elmgeb.f90
    entity.f90
    ftlay.f90
    getxyz.f90
    header.f90
    indxf.f90
    laytab.f90
    line.f90
    lwpolyline.f90
    matgeb.f90
    mtext.f90
    objects.f90
    outdxf.f90
    pline.f90
    putlin.f90
    putnod.f90
    rbtxt.f90
    readky.f90
    sectio.f90
    setzwg.f90
    snzwg.f90
    sortnt.f90
    tables.f90
    testac.f90
    text.f90
    unsect.f90
    vertex.f90
    )
endif()

# Additional files for the current project are stored in variable SRC2
set(SRC2 
  ${CMAKE_SOURCE_DIR}/Source/Post/circle.f90
  ${CMAKE_SOURCE_DIR}/Source/Post/struc.f90
  ${CMAKE_SOURCE_DIR}/Source/Post/zanfpp.f90
)

#  L I N U X
# Set the file in the variable for the linux version of acadini.f90
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(SRC3 ${CMAKE_SOURCE_DIR}/Source/System/Linux/acadnet/acadini_linux.f90)
endif()

# Add an executable to the project using the source files from SRC and SRC2 and/or acadnet
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_executable(Acadnet WIN32 ${SRC} ${SRC2})
else()
  add_executable(acadnet ${SRC} ${SRC2} ${SRC3})
endif()

# Add a dependency according to the OS
# On Linux, libpgplot.so is to be installed beforehand

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# Add a dependency between targets (target-name depend-target1 depend-target2 ...)
   add_dependencies(acadnet Femlib System)
# Link a target to given libraries 
   target_link_libraries(acadnet Femlib pgplot X11 System)
else()
   add_dependencies(Acadnet Types_Interfaces Femlib Femlib2D Pgplot System)
   target_link_libraries(Acadnet Femlib Pgplot System)
endif()

#  W I N D O W S
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# Intel Composer
    if(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel"))
# Additional compiler flags
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /libs:qwin")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /check:all /warn:all")
        set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} /check:bounds,pointers,uninit /warn:all")
# PGI Compiler
    elseif(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI"))
# Flang Compiler
    elseif(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang"))
# GFortran (GNU GCC) Compiler
    elseif(("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU"))
    endif() # Intel 
endif() # Windows


#  L I N U X
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set_target_properties(acadnet PROPERTIES COMPILE_FLAGS "-fpscomp filesfromcmd -warn all")
endif()

#Install the libraries and/or executables to the CMAKE_INSTALL_PRFIX path
install(TARGETS acadnet
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)