cmake_minimum_required(VERSION 2.8)

project("sycomore")

set(sycomore_MAJOR_VERSION 0)
set(sycomore_MINOR_VERSION 1)
set(sycomore_PATCH_VERSION 0)
set(sycomore_VERSION
    ${sycomore_MAJOR_VERSION}.${sycomore_MINOR_VERSION}.${sycomore_PATCH_VERSION})

option(BUILD_SHARED_LIBS "Build Sycomore with shared libraries." ON)
option(BUILD_TESTING "Build unit tests." ON)
option(BUILD_EXAMPLES "Build the examples directory." ON)

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

file(GLOB_RECURSE docs mkdocs.yml docs/*)
add_custom_target(
    Documentation ${CMAKE_COMMAND} -E echo "Documentation"
    SOURCES README.md ${docs})
set_target_properties(Documentation PROPERTIES FOLDER "Utils")

file(GLOB_RECURSE ContinuousIntegration .gitlab.yml .ci/*)
add_custom_target(
    ContinuousIntegration ${CMAKE_COMMAND} -E echo "Continuous Integration"
    SOURCES ${ContinuousIntegration})
set_target_properties(ContinuousIntegration PROPERTIES FOLDER "Utils")

add_subdirectory("src")

 if(BUILD_TESTING)
     add_subdirectory("tests")
 endif()

 if(BUILD_EXAMPLES)
     add_subdirectory("examples")
 endif()
