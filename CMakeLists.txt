cmake_minimum_required(VERSION 3.16)
project(QCProject LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build types
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release or Debug)" FORCE)
endif()

# Compiler flags per build type
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -fno-openmp")

# Set SIMD explicitly for Release (optional, already in -O3)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")

# OpenMP setup
find_package(OpenMP)
if(OpenMP_CXX_FOUND AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(${OpenMP_CXX_FLAGS})
    link_libraries(${OpenMP_CXX_LIBRARIES})
endif()

# BLAS/LAPACK backend selection
option(USE_MKL "Use Intel MKL" OFF)
option(USE_OPENBLAS "Use OpenBLAS" OFF)
option(USE_AMD_BLAS "Use AMD BLAS/libFLAME" OFF)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(common/include)

# Subdirectories
add_subdirectory(blasbind)
add_subdirectory(libmatrix)
add_subdirectory(libqcex)
add_subdirectory(qcint)