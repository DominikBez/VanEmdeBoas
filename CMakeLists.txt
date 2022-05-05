cmake_minimum_required (VERSION 2.8)
cmake_policy(SET CMP0048 NEW)
project(VanEmdeBoas)

#### USER DEFINED ##############################################################

#### BASIC SETTINGS ############################################################

include_directories(.)

set (CMAKE_CXX_FLAGS "-std=gnu++17 -pthread -march=native -Wall -Wextra -O3")
# -std=gnu++xx necessary for proper __int128 support (otherwise type traits do not work correctly)

#### TARGETS ###################################################################

set(GTEST_ROOT /usr/src/googletest)
find_path(GTEST_ROOT googletest/CMakeLists.txt HINTS ${GTEST_ROOT})

if (NOT ${GTEST_ROOT} STREQUAL GTEST_ROOT-NOTFOUND)
    enable_testing()

    add_subdirectory(${GTEST_ROOT} ${CMAKE_BINARY_DIR}/googletest EXCLUDE_FROM_ALL)

    add_executable(veb test/test_van_emde_boas.cpp)
    target_link_libraries(veb gtest_main)
	set_target_properties(veb PROPERTIES COMPILE_FLAGS "-g")
	
	add_executable(parallel_veb test/parallel_test_van_emde_boas.cpp)
	target_link_libraries(parallel_veb gtest_main)
	set_target_properties(parallel_veb PROPERTIES COMPILE_FLAGS "-g")

    include(GoogleTest)
    gtest_discover_tests(veb)
	gtest_discover_tests(parallel_veb)
endif()

add_executable(timing_veb test/timing_van_emde_boas.cpp)
set_target_properties(timing_veb PROPERTIES COMPILE_FLAGS "-DNDEBUG")

find_package(OpenMP REQUIRED)
add_executable(parallel_timing_veb test/parallel_timing_van_emde_boas.cpp)
set_target_properties(parallel_timing_veb PROPERTIES COMPILE_FLAGS "-DNDEBUG")
target_link_libraries(parallel_timing_veb PRIVATE OpenMP::OpenMP_CXX)
