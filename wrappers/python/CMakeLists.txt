if(NOT MSVC)
    set(XTENSOR_USE_XSIMD 1)
endif()

find_package(Python COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(pybind11 REQUIRED)
find_package(xsimd REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-python REQUIRED)

file(GLOB_RECURSE header_files "*.h")
file(GLOB_RECURSE template_files "*.txx")
file(GLOB_RECURSE source_files "*.cpp")
file(GLOB_RECURSE python_files "*.py")
list(SORT header_files)
list(SORT source_files)
list(SORT python_files)

pybind11_add_module(pysycomore SHARED ${source_files} ${header_files})

target_include_directories(
    pysycomore 
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/src ${PYTHON_INCLUDE_DIRS} ${xsimd_INCLUDE_DIRS})

target_link_libraries(
    pysycomore PUBLIC libsycomore Python::NumPy xtensor xtensor-python)

set_target_properties(
    pysycomore PROPERTIES 
    OUTPUT_NAME _sycomore
    $<$<PLATFORM_ID:Darwin>:SUFFIX .so>)

add_custom_target(
    pysycomore-pure ${CMAKE_COMMAND} -E echo "Pure-python files"
    SOURCES ${python_files})

execute_process(
    COMMAND ${Python_EXECUTABLE}
        -c "import os; import sysconfig; \
            scheme = {'posix': 'posix_prefix', 'nt': 'nt'}[os.name]; \
            path = sysconfig.get_path('purelib', scheme, {'base': '.'}); \
            print(path.replace(os.path.sep, '/'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE)

install(DIRECTORY DESTINATION "${PYTHON_SITE_PACKAGES}")
install(TARGETS pysycomore DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")
install(FILES ${python_files} DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")

if(WIN32)
    install(TARGETS libsycomore DESTINATION "${PYTHON_SITE_PACKAGES}/sycomore")
endif()

foreach(include_file ${header_files} ${template_files})
    file(
        RELATIVE_PATH 
        include_path ${CMAKE_SOURCE_DIR} ${include_file})
    get_filename_component(include_path ${include_path} PATH)
    install(FILES ${include_file} DESTINATION "include/sycomore/${include_path}")
endforeach()
