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
    
    if(MSVC)
        set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/src/sycomore/epg/simd_api_sse2.cpp PROPERTIES
            COMPILE_FLAGS "/arch:SSE2")
        set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/src/sycomore/epg/simd_api_avx.cpp PROPERTIES
            COMPILE_FLAGS "/arch:AVX")
        set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/src/sycomore/epg/simd_api_avx512.cpp PROPERTIES
            COMPILE_FLAGS "/arch:AVX512")
    else()
        set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/src/sycomore/epg/simd_api_sse2.cpp PROPERTIES
            COMPILE_FLAGS "-msse2")
        set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/src/sycomore/epg/simd_api_avx.cpp PROPERTIES
            COMPILE_FLAGS "-mavx")
        
        # gcc -march=skylake-avx512 -Q --help=target | grep avx512 | grep enabled
        # This is the earliest architecture to support AVX512
        set(AVX512_FLAGS "-mavx512bw -mavx512cd -mavx512dq -mavx512f -mavx512vl")
        if(APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            # Pretend we are GCC >= 6 so that xsimd enables AVX512
            # cf. https://github.com/xtensor-stack/xsimd/blob/7.4.9/include/xsimd/config/xsimd_instruction_set.hpp#L109
            set(AVX512_FLAGS "${AVX512_FLAGS} -U__GNUC__ -D__GNUC__=6")
        endif()
        set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/src/sycomore/epg/simd_api_avx512.cpp PROPERTIES 
            COMPILE_FLAGS "${AVX512_FLAGS}")
    endif()
    
    set(source_files ${source_files} ${lib_files})
endif()

pybind11_add_module(pysycomore SHARED ${source_files} ${header_files})

target_compile_definitions(
    pysycomore
    PRIVATE
    $<$<AND:$<PLATFORM_ID:Windows>,$<BOOL:${BUILD_STANDALONE_PYTHON_WRAPPERS}>>:BUILDING_SYCOMORE>)

target_include_directories(
    pysycomore 
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/src ${PYTHON_INCLUDE_DIRS} ${xsimd_INCLUDE_DIRS})

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    target_link_libraries(pysycomore PUBLIC libsycomore)
endif()

set_target_properties(
    pysycomore PROPERTIES 
    OUTPUT_NAME _sycomore
    $<$<PLATFORM_ID:Darwin>:SUFFIX .so>)

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
