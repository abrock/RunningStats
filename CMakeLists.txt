cmake_minimum_required(VERSION 2.6)
project(RunningStats)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -Werror -pedantic")

# Set variables
include(${CMAKE_SOURCE_DIR}/cmake/SetEnv.cmake)

# Library sources
add_subdirectory(${LIBRARY_FOLDER})

# Library examples
add_subdirectory(test)

# Install targets
include(${CMAKE_SOURCE_DIR}/cmake/InstallConfig.cmake)

add_executable(runningstats-autohist runningstats-autohist.cpp)

target_link_libraries(runningstats-autohist
    runningstats
    )

install(TARGETS runningstats-autohist DESTINATION bin)
