rock_library(sonaroctomap
    SOURCES SonarOcTree.cpp
    HEADERS SonarOcTree.hpp
    LIBS libmatio.so.2
    DEPS_CMAKE octomap)

rock_executable(sonaroctomap_bin Main.cpp
    LIBS glut GLU GL m  pthread
    DEPS sonaroctomap)

rock_executable(sonaroctomap_bin_extra MainExtra.cpp
    DEPS sonaroctomap)

