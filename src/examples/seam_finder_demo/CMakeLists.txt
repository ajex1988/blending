#
# Copyright (c) 2018 loyave
#

# Target name
set(target seam_finder_demo)

# Find Package
find_package(OPENCV REQUIRED)

# Includes
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

# Sources
file(GLOB sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

# Build executable
add_executable(${target}
    ${sources})

# Project options
set_target_properties(${target}
    PROPERTIES
    ${DEFAULT_PROJECT_OPTIONS}
)

# Compile options
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_definitions(-D_USE_MATH_DEFINES) # for M_PI
endif ()

target_compile_options(${target}
    PRIVATE

    PUBLIC
    ${DEFAULT_COMPILE_OPTIONS}

    INTERFACE
)

target_compile_definitions(${target}
    PRIVATE
    RESOURCES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../../resources/"
)

# Link libraries
target_link_libraries(${target}
    PRIVATE
    ${DEFAULT_LINKER_OPTIONS}
    ${OpenCV_LIBS}
    utils
    blending)
