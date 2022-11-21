# Subdirectories in which other CMakeLists.txt or code files are located
# For Linux use Blas in Intel's MKL

add_subdirectory(Geompack)
add_subdirectory(Json)
add_subdirectory(Umfpack)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_subdirectory(Blas)
    add_subdirectory(Sparsekit)
    add_subdirectory(Pgplot)
    add_subdirectory(Umfpack55)
endif()