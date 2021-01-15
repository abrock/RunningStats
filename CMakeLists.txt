cmake_minimum_required(VERSION 2.6)
project(RunningStats)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -Werror -pedantic")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set variables
include(${CMAKE_SOURCE_DIR}/cmake/SetEnv.cmake)

# Library sources
add_subdirectory(${LIBRARY_FOLDER})

# Library examples
add_subdirectory(test)

# Install targets
include(${CMAKE_SOURCE_DIR}/cmake/InstallConfig.cmake)

foreach(PROGRAM runningstats-autohist runningstats-print)
add_executable(${PROGRAM} ${PROGRAM}.cpp)

target_link_libraries(${PROGRAM}
    runningstats
    )

install(TARGETS ${PROGRAM} DESTINATION bin)
endforeach()
