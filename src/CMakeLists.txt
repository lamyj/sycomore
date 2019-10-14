find_package(PythonInterp REQUIRED)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    file(GLOB_RECURSE header_files "*.h")
    file(GLOB_RECURSE source_files "*.cpp")
    file(GLOB_RECURSE template_files "*.txx")
    list(SORT header_files)
    list(SORT source_files)
    list(SORT template_files)

    include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
    
    if(WIN32)
        add_definitions(-DBUILDING_SYCOMORE)
    endif()

    add_library(libsycomore ${source_files} ${header_files} ${template_files})
    if(APPLE AND OPENMP_FOUND)
        target_link_libraries(libsycomore omp)
    endif()
    set_target_properties(libsycomore PROPERTIES OUTPUT_NAME sycomore)
    set_target_properties(libsycomore PROPERTIES
        VERSION ${sycomore_VERSION}
        SOVERSION ${sycomore_MAJOR_VERSION})

    install(
        DIRECTORY sycomore/ DESTINATION include/sycomore
        FILES_MATCHING PATTERN "*.h" PATTERN "*.txx")
    install(
        TARGETS libsycomore
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
endif()

file(GLOB_RECURSE python_files "*.py")
list(SORT python_files)

add_custom_target(
    libsycomore-python ${CMAKE_COMMAND} -E echo "Python files"
    SOURCES ${python_files})

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE}
          -c "from distutils.sysconfig import *; print(get_python_lib(True, prefix=''))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
    set(PYTHON_SITE_PACKAGES ".")
endif()

foreach(python_file ${python_files})
    file(RELATIVE_PATH module ${CMAKE_CURRENT_SOURCE_DIR} ${python_file})
    get_filename_component(module ${module} PATH)
    install(FILES ${python_file} DESTINATION "${PYTHON_SITE_PACKAGES}/${module}")
endforeach()
