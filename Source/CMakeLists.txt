#Subdirectories in which other CMakeLists.txt or code files are located

add_subdirectory(Femlib)
add_subdirectory(Femlib3D)
add_subdirectory(Post3D)
add_subdirectory(Solver3D)
add_subdirectory(System)

#  W I N D O W S
# Only include in the windows version
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    IF(("${CMAKE_GENERATOR}" MATCHES "Visual Studio"))
        add_subdirectory(Acadnet)
        add_subdirectory(Triang)
        add_subdirectory(Fluidlib)
        add_subdirectory(Solver)
        add_subdirectory(Post)
        add_subdirectory(Userinterface)
        add_subdirectory(Userinterface3D)
    ENDIF()
ENDIF()


# L I N U X
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_subdirectory(Acadnet)
  add_subdirectory(Triang)
  add_subdirectory(Fluidlib)
  add_subdirectory(Solver)
  add_subdirectory(Post)
  add_subdirectory(Userinterface)
  add_subdirectory(Userinterface3D)
endif()