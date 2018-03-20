#
# Copyright (c) 2018 loyave
#

# CMake version
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

# Include cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Declare project
project(blending)

# Set output directories
set(DEFAULT_CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Includes
include_directories(include)

# Compile options
include(cmake/CompileOptions.cmake)

# Build type - Release by default
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Unpack resources
set(RESOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/resources)
set(RESOURCES_OBJS
        ${RESOURCES_DIR}/script.txt
        ${RESOURCES_DIR}/0.png
        ${RESOURCES_DIR}/1.png
        ${RESOURCES_DIR}/2.png
        ${RESOURCES_DIR}/3.png
        ${RESOURCES_DIR}/4.png
        ${RESOURCES_DIR}/5.png)

add_subdirectory(src/blending)
add_subdirectory(src/utils)

add_subdirectory(src/examples/remapper_demo)
add_subdirectory(src/examples/seam_finder_demo)
add_subdirectory(src/examples/mvc_blender_demo)
add_subdirectory(src/examples/multiband_blender_demo)
#add_subdirectory(src/examples/fb_blend_demo)