find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

rock_library(sonaroctomap
    SOURCES SonarOcTree.cpp
    HEADERS SonarOcTree.hpp
    DEPS_PKGCONFIG matio
    DEPS_CMAKE octomap)

rock_executable(sonaroctomap_bin Main.cpp
    LIBS pthread
    DEPS sonaroctomap)

rock_executable(sonaroctomap_bin_extra MainExtra.cpp
    DEPS sonaroctomap)

