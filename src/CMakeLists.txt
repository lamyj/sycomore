find_package(Eigen3 REQUIRED)

file(GLOB_RECURSE header_files "*.h")
file(GLOB_RECURSE source_files "*.cpp")
file(GLOB_RECURSE template_files "*.txx")
list(SORT header_files)
list(SORT source_files)
list(SORT template_files)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(libsycomore ${source_files} ${header_files} ${template_files})
target_link_libraries(libsycomore profiler)
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
