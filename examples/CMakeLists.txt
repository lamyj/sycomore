option(INSTALL_EXAMPLES "Install the examples" OFF)

if(NOT BUILD_STANDALONE_PYTHON_WRAPPERS)
    add_subdirectory("cpp")
endif()
