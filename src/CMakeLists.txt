rock_library(sonaroctomap
    SOURCES SonarOcTree.cpp
    HEADERS SonarOcTree.hpp
    LIBS libmatio.so
    DEPS_CMAKE octomap)

rock_executable(sonaroctomap_bin Main.cpp
    DEPS sonaroctomap)

