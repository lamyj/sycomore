if(NOT MSVC)
    set(XTENSOR_USE_XSIMD 1)
endif()

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
find_package(xsimd REQUIRED)
find_package(xtensor REQUIRED)

include_directories(
    ${CMAKE_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS} ${xsimd_INCLUDE_DIRS})

if(BUILD_SHARED_LIBS)
    add_definitions(-D BOOST_ALL_DYN_LINK)
endif()

file(GLOB_RECURSE cpp_tests cpp/*.cpp)

foreach(test_file ${cpp_tests})
    get_filename_component(unit ${test_file} NAME_WE)
    get_filename_component(module ${test_file} PATH)
    get_filename_component(module ${module} NAME)
    if(NOT "${module}" STREQUAL "cpp")
        set(test "${module}_${unit}")
    else()
        set(test "${unit}")
    endif()

    add_executable(test_${test} ${test_file})
    
    target_compile_definitions(
        test_${test} PRIVATE $<$<BOOL:BUILD_SHARED_LIBS>:BOOST_ALL_DYN_LINK>)
    
    target_include_directories(
        test_${test} PRIVATE ${CMAKE_SOURCE_DIR}/src ${Boost_LIBRARY_DIRS})
    
    target_link_libraries(
        test_${test} PRIVATE ${Boost_LIBRARIES} libsycomore xtensor)
    
    set_target_properties(test_${test} PROPERTIES OUTPUT_NAME ${test})

    file(READ ${test_file} content)

    set(pattern "BOOST_(AUTO|FIXTURE)_TEST_CASE\\(([^),]+)")

    string(REGEX MATCHALL ${pattern} cases ${content})

    foreach(case ${cases})
        string(REGEX REPLACE ${pattern} "\\2" case ${case})
        if(NOT ${case} MATCHES "##")
            add_test("${test}_${case}" "${test}" "--run_test=${case}")
        endif()
    endforeach()
endforeach()

file(GLOB_RECURSE data data/*)
add_custom_target(data ${CMAKE_COMMAND} -E echo "Data" SOURCES ${data})
set_target_properties(data PROPERTIES FOLDER "Test data")

file(GLOB_RECURSE python_tests python/*.py)
add_custom_target(
    python_tests ${CMAKE_COMMAND} -E echo "Python tests" SOURCES ${python_tests})
