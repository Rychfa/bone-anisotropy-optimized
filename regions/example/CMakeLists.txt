cmake_minimum_required(VERSION 3.13)
project(example)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
#set(CMAKE_VERBOSE_MAKEFILE "ON")

include_directories(
    ../
    ../../utils
)

set(SOURCE_FILES
        simple_memory.c
        )

set(CMAKE_C_FLAGS "-O3 -mavx2 -mfma -fno-tree-vectorize")

add_executable(regions_extraction ${SOURCE_FILES})

# Add math library for floating point environment
#target_link_libraries(example m)