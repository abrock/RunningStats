cmake_minimum_required(VERSION 2.6)
project(RunningStats)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector-strong -Wall -Werror=shadow -pedantic")

OPTION(PRODUCTION "Enable this for production builds with less strict compile error checking" ON)
if (PRODUCTION)
else()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(fmt REQUIRED)
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(TCLAP REQUIRED)

# Set variables
include(${CMAKE_SOURCE_DIR}/cmake/SetEnv.cmake)

# Library sources
add_subdirectory(${LIBRARY_FOLDER})

option(BUILD_TESTS "" OFF)
if(BUILD_TESTS)
# Library examples
add_subdirectory(test)
endif()

# Install targets
include(${CMAKE_SOURCE_DIR}/cmake/InstallConfig.cmake)

foreach(PROGRAM runningstats-autohist runningstats-print runningstats-thresholds runningstats-get-median)
add_executable(${PROGRAM} ${PROGRAM}.cpp)

target_link_libraries(${PROGRAM}
    runningstats
    )

install(TARGETS ${PROGRAM} DESTINATION bin)
endforeach()

option(USE_OPENCV "Use OpenCV to build an example program which reads images and plots pixel values")
if (USE_OPENCV)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(runningstats-plot-image plot-image.cpp)
target_link_libraries(runningstats-plot-image
    runningstats
    ${OpenCV_LIBRARIES}
    )
install(TARGETS runningstats-plot-image DESTINATION bin)
endif(USE_OPENCV)
