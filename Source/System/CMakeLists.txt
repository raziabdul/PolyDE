#Goto subdirectory PC if using Windows
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_subdirectory(PC)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#Goto subdirectory Linux if using Linux 
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_subdirectory(Linux)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")