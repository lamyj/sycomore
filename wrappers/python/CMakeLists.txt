find_package(pybind11 REQUIRED)

file(GLOB_RECURSE header_files "*.h")
file(GLOB_RECURSE source_files "*.cpp")
file(GLOB_RECURSE python_files "*.py")
list(SORT header_files)
list(SORT source_files)
list(SORT python_files)

include_directories(${CMAKE_SOURCE_DIR}/src ${PYTHON_INCLUDE_DIRS})

pybind11_add_module(pysycomore SHARED ${source_files} ${header_files})

set_target_properties(pysycomore PROPERTIES OUTPUT_NAME _sycomore)
if(APPLE)
    set_target_properties(pysycomore PROPERTIES SUFFIX ".so")
endif()

target_link_libraries(pysycomore ${JsonCpp_LIBRARIES} libsycomore)

add_custom_target(
    pysycomore-pure ${CMAKE_COMMAND} -E echo "Pure-python files"
    SOURCES ${python_files})


execute_process(
    COMMAND ${PYTHON_EXECUTABLE}
      -c "from distutils.sysconfig import *; print(get_python_lib(True, prefix=''))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE)

install(DIRECTORY DESTINATION "${PYTHON_SITE_PACKAGES}")
install(TARGETS pysycomore DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")
install(FILES ${python_files} DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")
