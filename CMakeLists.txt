cmake_minimum_required(VERSION 2.8)

project("sycomore")

set(sycomore_VERSION 1.0.0)
string(REGEX MATCH "^[0-9]+" sycomore_MAJOR_VERSION "${sycomore_VERSION}")

option(BUILD_SHARED_LIBS "Build Sycomore with shared libraries." ON)
option(BUILD_TESTING "Build unit tests." ON)
option(BUILD_PYTHON_WRAPPERS "Build the Python Wrappers." ON)
option(
    BUILD_STANDALONE_PYTHON_WRAPPERS
    "Build the Python Wrappers as a standalone library." OFF)
option(USE_OPENMP "Compile with OpenMP support." ON)

if(NOT BUILD_PYTHON_WRAPPERS AND BUILD_STANDALONE_PYTHON_WRAPPERS)
    message(WARNING "Standalone Python wrappers requested: setting BUILD_PYTHON_WRAPPERS to ON")
    set(BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build the Python Wrappers." FORCE)
endif()

if(USE_OPENMP AND BUILD_STANDALONE_PYTHON_WRAPPERS)
    message(WARNING "Standalone Python wrappers requested: setting USE_OPENMP to OFF")
    set(USE_OPENMP OFF CACHE BOOL "Compile with OpenMP support." FORCE)
endif()

include(CTest)

# Add the C++0x or C++11 flag
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG(-std=c++0x COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()

if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif()

file(GLOB_RECURSE docs docs/*)
add_custom_target(
    Documentation ${CMAKE_COMMAND} -E echo "Documentation"
    SOURCES README.md ${docs})

file(GLOB_RECURSE ContinuousIntegration .gitlab-ci.yml .ci/*)
add_custom_target(
    ContinuousIntegration ${CMAKE_COMMAND} -E echo "Continuous Integration"
    SOURCES ${ContinuousIntegration})
set_target_properties(ContinuousIntegration PROPERTIES FOLDER "Utils")

add_subdirectory("src")

if(BUILD_PYTHON_WRAPPERS)
    add_subdirectory("wrappers/python")
endif()

if(BUILD_TESTING)
    add_subdirectory("tests")
endif()
