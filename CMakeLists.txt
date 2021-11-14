cmake_minimum_required(VERSION 3.5)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(SET CMP0074 NEW)
endif()

project("sycomore" VERSION 1.3.2)

option(BUILD_SHARED_LIBS "Build Sycomore with shared libraries." ON)
option(BUILD_TESTING "Build unit tests." ON)
option(BUILD_PYTHON_WRAPPERS "Build the Python Wrappers." ON)
option(
    BUILD_STANDALONE_PYTHON_WRAPPERS
    "Build the Python Wrappers as a standalone library." OFF)

if(NOT BUILD_PYTHON_WRAPPERS AND BUILD_STANDALONE_PYTHON_WRAPPERS)
    message(WARNING "Standalone Python wrappers requested: setting BUILD_PYTHON_WRAPPERS to ON")
    set(BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build the Python Wrappers." FORCE)
endif()

set(CMAKE_INSTALL_MESSAGE LAZY)

include(CTest)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# Export the build tree (don't install the generated file)
if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    export(
        TARGETS libsycomore NAMESPACE sycomore:: 
        FILE "${PROJECT_BINARY_DIR}/sycomoreTargets.cmake")

    # Export the install tree
    write_basic_package_version_file(
        "sycomoreConfigVersion.cmake" COMPATIBILITY SameMajorVersion)
    set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}/lib/CMake/sycomore")
    install(EXPORT sycomoreTargets NAMESPACE sycomore:: DESTINATION "${INSTALL_CMAKE_DIR}")
    configure_file(sycomoreConfig.cmake.in sycomoreConfig.cmake @ONLY)
    install(
        FILES 
            ${CMAKE_CURRENT_BINARY_DIR}/sycomoreConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/sycomoreConfigVersion.cmake
        DESTINATION "${INSTALL_CMAKE_DIR}")
endif()
