cmake_minimum_required(VERSION 3.0)

if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

#option(GLPK_WITH_GMP "use GNU MP bignum library" OFF)
set(GLPK_WITH_GMP ON)

# Define variables for library
set(PROJECT_NAME glpk-simplex-trace)
set(LIBRARY_NAME glpk-simplex-trace)

project(${PROJECT_NAME} C)

enable_testing()
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type." FORCE)
endif ()

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/src)

# add_subdirectory("examples")
add_subdirectory("src")
