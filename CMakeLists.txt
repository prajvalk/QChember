cmake_minimum_required(VERSION 3.10)
project(QChember LANGUAGES C CXX)

# Allow user to set build type if not set
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Defaulting to RelWithDebInfo.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
endif()

# Expose our helper modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CompilerHelpers)
include(FindCUDA_CUBLAS)
include(FindBLAS_LAPACK)

# Set sane default output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Configure default flags for the detected compiler
configure_compiler_flags()

# Detect CUDA / CUBLAS
#detect_cuda_and_cublas()
detect_blas_lapack()

include_directories(include)

add_library(qchember SHARED modules/lacis/lacis.cpp
                            modules/lacis/matrixmarketio.cpp)

# Example: create a library target and attach common flags
#add_library(example_lib src/example.cpp)
#target_link_libraries(example_lib PRIVATE ${CMAKE_THREAD_LIBS_INIT})
#apply_common_target_flags(example_lib)

# If CUDA exists, create an interface target to link CUDA libs
#if(HAS_CUDA)
 # message(STATUS "CUDA detected: enabling CUDA-related targets")
 # add_library(project_cuda INTERFACE)
 # if(HAS_CUBLAS)
 #   target_link_libraries(project_cuda INTERFACE ${CUBLAS_LIBRARIES})
 # endif()
 # target_compile_definitions(project_cuda INTERFACE HAS_CUDA=1)
#endif()