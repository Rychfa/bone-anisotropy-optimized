cmake_minimum_required(VERSION 3.13)
project(bonemap)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

# ================================================================================
# Reset the C/C++ flags, and take full control
# ================================================================================
set(CMAKE_C_FLAGS           "")
set(CMAKE_CXX_FLAGS         "")
set(CMAKE_C_FLAGS_RELEASE   "")
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_C_FLAGS_DEBUG     "")
set(CMAKE_CXX_FLAGS_DEBUG   "")

# ================================================================================
# Include directories and set source files
# ================================================================================
include_directories(
        reader
        mil
        regions
        ellipsoid
        utils
        eigen
)

set(SOURCE_FILES
        main.cpp
        mil/mil.c
        mil/test_mil.c
        ellipsoid/ellipsoid.c
        reader/reader.c
        regions/region_extraction.c
        eigen/eigen.c
        utils/utils.c
        kernel.c)


# ================================================================================
# C++ Global Flags
# ================================================================================
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

    # clang - C++ flags
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    # GNU GCC - C++ flags
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

    # Intel ICC - C++ flags
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

    # Microsoft Visual Studio - C++ flags
    set(CMAKE_CXX_FLAGS "-O2 -std=c++11")

endif()

# ================================================================================
# C Global Flags
# ================================================================================

if (CMAKE_C_COMPILER_ID MATCHES "Clang")

    # clang - CC flags
    set(CMAKE_C_FLAGS "-O3")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_C_FLAGS "-O3")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

    # Intel ICC - CC flags
    set(CMAKE_C_FLAGS "-O3 -xHost")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")

    # Microsoft Visual Studio - CC flags
    set(CMAKE_C_FLAGS "-O2 -std=c++11")

endif()

# ================================================================================
# Define binary executable and libraries to be linked
# ================================================================================
add_executable(bonemap ${SOURCE_FILES})

# Add math library for floating point environment
target_link_libraries(bonemap m)