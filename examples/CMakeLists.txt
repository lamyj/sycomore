file(GLOB_RECURSE examples *.cpp)

foreach(example_file ${examples})
    get_filename_component(example ${example_file} NAME_WE)
    add_executable(${example} ${example_file})
    target_link_libraries(${example} libsycomore)
endforeach()
