find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

file(GLOB_RECURSE header_files "*.h")
file(GLOB_RECURSE source_files "*.cpp")
file(GLOB_RECURSE template_files "*.txx")
list(SORT header_files)
list(SORT source_files)
list(SORT template_files)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})

add_library(libsycomore ${source_files} ${header_files} ${template_files})
if(APPLE)
    target_link_libraries(libsycomore omp)
endif()
set_target_properties(libsycomore PROPERTIES OUTPUT_NAME sycomore)
set_target_properties(libsycomore PROPERTIES
    VERSION ${sycomore_VERSION}
    SOVERSION ${sycomore_MAJOR_VERSION})

install(
    DIRECTORY sycomore/ DESTINATION include/sycomore
    FILES_MATCHING PATTERN "*.h")
install(
    TARGETS libsycomore
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
