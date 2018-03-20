#
# Copyright (c) 2018 loyave
#

# Target name
set(target utils)

# Define
set(root_dir ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Find Package
find_package(OPENCV REQUIRED)
find_package(CGAL REQUIRED)

# Includes
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)


# Sources
file(GLOB header_dir
    ${root_dir}/include/${target})

file(GLOB headers
    ${header_dir}/*.h)

file(GLOB sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

# Build library
add_library(${target}
    ${sources})

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
