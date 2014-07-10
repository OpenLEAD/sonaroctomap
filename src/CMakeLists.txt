rock_library(sonaroctomap
    SOURCES SonarOcTree.cpp
    HEADERS SonarOcTree.h 
    DEPS_CMAKE octomap)

rock_executable(sonaroctomap_bin Main.cpp
    DEPS sonaroctomap)

