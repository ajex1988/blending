#
# Copyright (c) 2018 loyave
#

# Target name
set(target blending)

# Define
set(root_dir ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Find Package
find_package(OPENCV REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)

# Includes

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)
include(${CGAL_USE_FILE})

# Sources
file(GLOB header_dir
    ${root_dir}/include/${target})

file(GLOB headers
    ${header_dir}/*.h)

file(GLOB sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

# Custom-build event
set(blending_header_gen_py ${root_dir}/scripts/header_gen.py)
add_custom_command(
    COMMAND python ${blending_header_gen_py}
    DEPENDS ${headers}
    OUTPUT ${root_dir}/include/blending.h
)
add_custom_target(blending_header_gen_py ALL
    DEPENDS ${root_dir}/include/blending.h)

# Build library
add_library(${target}
    ${sources})
add_dependencies(${target} blending_header_gen_py)

# Project options
set_target_properties(${target}
    PROPERTIES
    ${DEFAULT_PROJECT_OPTIONS}
)

# Compile options
target_compile_options(${target}
    PRIVATE

    PUBLIC
    ${DEFAULT_COMPILE_OPTIONS}

    INTERFACE
)

target_link_libraries(${target}
    PRIVATE

    PUBLIC
    ${DEFAULT_LINKER_OPTIONS}
    ${DEFAULT_LIBRARIES}    
    ${OpenCV_LIBS}
    INTERFACE
)

# Install
install(TARGETS ${target} DESTINATION lib)
install(DIRECTORY ${header_dir} DESTINATION include)
