cmake_minimum_required(VERSION 3.31)
project(newscf_cis)

set(CMAKE_CXX_STANDARD 23)


# ================================
# Options
# ================================
option(USE_OPENMP  "Enable OpenMP parallelization" ON)
option(USE_LIBCINT "Build with external libcint module" ON)

# ================================
# Set Build Types
# ================================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# ================================
# Compiler Flags
# ================================
set(COMMON_CXX_FLAGS "-Wall -Wextra -Wpedantic")
set(COMMON_C_FLAGS   "-Wall -Wextra -Wpedantic")

# Detect compiler
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS_DEBUG   "${COMMON_C_FLAGS} -g")
    set(CMAKE_CXX_FLAGS_DEBUG "${COMMON_CXX_FLAGS} -g")
    set(CMAKE_C_FLAGS_RELEASE   "${COMMON_C_FLAGS} -O3 -DNDEBUG -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_CXX_FLAGS} -O3 -DNDEBUG -march=native")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS_DEBUG   "${COMMON_C_FLAGS} -g")
    set(CMAKE_CXX_FLAGS_DEBUG "${COMMON_CXX_FLAGS} -g")
    set(CMAKE_C_FLAGS_RELEASE   "${COMMON_C_FLAGS} -O3 -DNDEBUG -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_CXX_FLAGS} -O3 -DNDEBUG -march=native")

elseif(MSVC)
    message(WARNING "This boilerplate currently doesn't optimize MSVC flags.")
endif()

# ================================
# OpenMP
# ================================
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message(STATUS "Using OpenMP")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

include_directories(public_api)
include_directories(private_api)

# Test Data
set(TEST_DATA_SOURCE_DIR "${CMAKE_SOURCE_DIR}/TestData")
set(TEST_DATA_DEST_DIR "${CMAKE_BINARY_DIR}/TestData")

add_custom_target(copy_test_data ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${TEST_DATA_SOURCE_DIR}
        ${TEST_DATA_DEST_DIR}
        COMMENT "Copying test data to ${TEST_DATA_DEST_DIR}"
)

enable_testing()

add_subdirectory(src/lacis)
add_subdirectory(src/lacislapack)
add_subdirectory(src/newscf-utils)

if(USE_LIBCINT)
    add_subdirectory(src/libcint-wrapper)
endif()

add_subdirectory(src/scf)