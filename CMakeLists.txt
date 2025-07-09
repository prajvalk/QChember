cmake_minimum_required(VERSION 3.16)
project(newscf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build types
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release or Debug)" FORCE)
endif()

# Compiler flags per build type
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -fno-openmp -DNEWSCF_LOG_LEVEL=DEV_DUMP")

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

include_directories(include)
enable_testing()

add_subdirectory(qcint)

file(GLOB_RECURSE NEWSCF_SRC src/*.cpp)
file(GLOB_RECURSE NEWSCF_TESTS tests/*.cpp)

add_library(newscf SHARED ${NEWSCF_SRC})
target_include_directories(newscf PUBLIC qcint/include)
target_link_libraries(newscf PUBLIC cint ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

file(GLOB_RECURSE NEWSCF_DATA data/*)
file(COPY ${NEWSCF_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

foreach (test_file ${NEWSCF_TESTS})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})

    target_include_directories(${test_name} PRIVATE include)
    target_link_libraries(${test_name} PRIVATE newscf cint ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${test_name} PRIVATE -Og -g3 -fno-openmp)
    else()
        target_compile_options(${test_name} PRIVATE -O3 -march=native)
        if(OpenMP_CXX_FOUND)
            target_compile_options(${test_name} PRIVATE ${OpenMP_CXX_FLAGS})
            target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endif()

    get_filename_component(test_name ${test_file} NAME_WE)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach ()