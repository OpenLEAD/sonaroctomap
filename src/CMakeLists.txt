rock_library(sonaroctomap
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(sonaroctomap_bin Main.cpp
    DEPS sonaroctomap)

