cmake_minimum_required(VERSION 3.0)

# Define variables for library
set(PROJECT_NAME glpk-debug)
set(LIBRARY_NAME glpk-debug)

project(${PROJECT_NAME} C)

enable_testing()
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type." FORCE)
endif ()

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/src)

# add_subdirectory("examples")
add_subdirectory("src")
