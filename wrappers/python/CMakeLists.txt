find_package(pybind11 REQUIRED)
find_package(xsimd REQUIRED)

file(GLOB_RECURSE header_files "*.h")
file(GLOB_RECURSE source_files "*.cpp")
file(GLOB_RECURSE python_files "*.py")
list(SORT header_files)
list(SORT source_files)
list(SORT python_files)

if(BUILD_STANDALONE_PYTHON_WRAPPERS)
    file(GLOB_RECURSE lib_files "${CMAKE_SOURCE_DIR}/src/*.cpp")
    set(source_files ${source_files} ${lib_files})
    if(WIN32)
        add_definitions(-DBUILDING_SYCOMORE)
    endif()
endif()

include_directories(
    ${CMAKE_SOURCE_DIR}/src ${PYTHON_INCLUDE_DIRS} ${xsimd_INCLUDE_DIRS})

pybind11_add_module(pysycomore SHARED ${source_files} ${header_files})

set_target_properties(pysycomore PROPERTIES OUTPUT_NAME _sycomore)
if(APPLE)
    set_target_properties(pysycomore PROPERTIES SUFFIX ".so")
endif()

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    target_link_libraries(pysycomore PUBLIC libsycomore)
endif()

add_custom_target(
    pysycomore-pure ${CMAKE_COMMAND} -E echo "Pure-python files"
    SOURCES ${python_files})

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE}
          -c "import os; \
            from distutils.sysconfig import *; \
            print(get_python_lib(True, prefix='').replace(os.path.sep, '/'))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
    set(PYTHON_SITE_PACKAGES ".")
endif()

install(DIRECTORY DESTINATION "${PYTHON_SITE_PACKAGES}")
install(TARGETS pysycomore DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")
install(FILES ${python_files} DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")
