find_package(PythonInterp REQUIRED)
find_package(xsimd REQUIRED)

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    file(GLOB_RECURSE header_files "*.h")
    file(GLOB_RECURSE source_files "*.cpp")
    file(GLOB_RECURSE template_files "*.txx")
    list(SORT header_files)
    list(SORT source_files)
    list(SORT template_files)
    
    if(MSVC)
        set_source_files_properties(
            sycomore/epg/simd_api_sse2.cpp PROPERTIES 
            COMPILE_FLAGS "/arch:SSE2")
        set_source_files_properties(
            sycomore/epg/simd_api_avx.cpp PROPERTIES
            COMPILE_FLAGS "/arch:AVX")
        set_source_files_properties(
            sycomore/epg/simd_api_avx512.cpp PROPERTIES
            COMPILE_FLAGS "/arch:AVX512")
    else()
        set_source_files_properties(
            sycomore/epg/simd_api_sse2.cpp PROPERTIES COMPILE_FLAGS "-msse2")
        set_source_files_properties(
            sycomore/epg/simd_api_avx.cpp PROPERTIES COMPILE_FLAGS "-mavx")
        
        # gcc -march=skylake-avx512 -Q --help=target | grep avx512 | grep enabled
        # This is the earliest architecture to support AVX512
        set(AVX512_FLAGS "-mavx512bw -mavx512cd -mavx512dq -mavx512f -mavx512vl")
        if(APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            # Pretend we are GCC >= 6 so that xsimd enables AVX512
            # cf. https://github.com/xtensor-stack/xsimd/blob/7.4.9/include/xsimd/config/xsimd_instruction_set.hpp#L109
            set(AVX512_FLAGS "${AVX512_FLAGS} -U__GNUC__ -D__GNUC__=6")
        endif()
        set_source_files_properties(
            sycomore/epg/simd_api_avx512.cpp PROPERTIES 
            COMPILE_FLAGS "${AVX512_FLAGS}")
    endif()
    
    add_library(libsycomore ${source_files} ${header_files} ${template_files})
    
    target_compile_definitions(
        libsycomore
        PRIVATE
            $<$<PLATFORM_ID:Windows>:BUILDING_SYCOMORE>)
    
    target_include_directories(
        libsycomore
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/> $<INSTALL_INTERFACE:>
            ${xsimd_INCLUDE_DIRS})
    
    set_target_properties(
        libsycomore PROPERTIES 
        OUTPUT_NAME sycomore 
        VERSION ${sycomore_VERSION} 
        SOVERSION ${sycomore_VERSION_MAJOR}
        $<$<PLATFORM_ID:Darwin>:INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib>)
    
    install(
        TARGETS libsycomore
        EXPORT sycomoreTargets
        ARCHIVE DESTINATION bin LIBRARY DESTINATION lib RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include/sycomore
        INCLUDES DESTINATION include)
    
    foreach(include_file ${header_files} ${template_files})
        file(
            RELATIVE_PATH 
            include_path ${CMAKE_CURRENT_SOURCE_DIR} ${include_file})
        get_filename_component(include_path ${include_path} PATH)
        install(FILES ${include_file} DESTINATION "include/${include_path}")
    endforeach()
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
