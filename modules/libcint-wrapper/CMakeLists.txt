project(newscf_libcint_wrapper)

include(FetchContent)

FetchContent_Declare(
        libcint
            GIT_REPOSITORY https://github.com/sunqm/qcint
        GIT_TAG v6.1.2
        GIT_SHALLOW TRUE
        GIT_PROGRESS ON
        SYSTEM
)

FetchContent_MakeAvailable(libcint)

add_library(newscf_libcint SHARED cint_init.cpp cint_hf.cpp)
target_link_libraries(newscf_libcint cint)
target_include_directories(newscf_libcint PRIVATE ${libcint_BINARY_DIR}/include)
target_include_directories(newscf_libcint PRIVATE ${libcint_SOURCE_DIR}/include)
target_include_directories(newscf_libcint PRIVATE ${libcint_BINARY_DIR}/src)
target_compile_options(newscf_libcint PRIVATE -march=native)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)