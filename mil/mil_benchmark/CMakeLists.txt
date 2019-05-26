#cmake_minimum_required(VERSION 3.13)
#project(mil_benchmark)
#set(CMAKE_C_STANDARD 11)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
##set(CMAKE_VERBOSE_MAKEFILE "ON")
#
#include_directories(
#    ../
#    ../../ellipsoid
#)
#
#set(SOURCE_FILES
#        main.cpp
#        ../mil2.c)
#
#add_executable(mil_benchmark ${SOURCE_FILES})
#
## Add math library for floating point environment
#target_link_libraries(mil_benchmark m)

##################################################################################

cmake_minimum_required(VERSION 3.13)
project(mil_benchmark)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_VERBOSE_MAKEFILE "ON")

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
../
../../ellipsoid
)

set(SOURCE_FILES
        main.cpp
        ../mil2.c
        ../mil2_simd.c
        mil_tests.c)


# ================================================================================
# C++ Global Flags
# ================================================================================
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

    # clang - C++ flags
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    # GNU GCC - C++ flags
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -fno-tree-vectorize -mavx2 -mfma")

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

    set(CMAKE_C_FLAGS "-O3 -fno-tree-vectorize -mavx2 -mfma")

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
add_executable(mil_benchmark ${SOURCE_FILES})

# Add math library for floating point environment
target_link_libraries(mil_benchmark m)