# (uncomment these lines for an external project)
#
# cmake_minimum_required(VERSION 2.6)
# project(Bar)
# find_package(Foo REQUIRED)
# include_directories(${FOO_INCLUDE_DIRS})


find_package(GTest)
if (GTest_FOUND)
find_package(Threads REQUIRED)

foreach(executable test-2d test-knuth)
add_executable(${executable} ${executable}.cpp)
endforeach()

add_executable(test-runningstats
    test-runningstats.cpp
    subtestimage.cpp
    subtestimage2d.cpp
    subtestcopyassign.cpp
    subteststatsnd.cpp
    subtest_weighted_runningstats.cpp
    subtest_threshold_error_mean.cpp
    subtest_precision.cpp
    subtest_quantilestats.cpp
    )

foreach(executable test-runningstats test-2d test-knuth)
target_link_libraries(${executable} ${LIBRARY_NAME})
target_link_libraries(${executable}
	GTest::gtest
	)
endforeach()

endif(GTest_FOUND)

#install(TARGETS test-runningstats
# In order to export target, uncomment next line
#   EXPORT ${PROJECT_EXPORT}
#  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
