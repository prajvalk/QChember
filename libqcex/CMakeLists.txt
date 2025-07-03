project(libqcex LANGUAGES CXX)

# Collect source files
file(GLOB_RECURSE LIBQCEX_SRC src/*.cpp)

# Create shared library
add_library(libqcex SHARED ${LIBQCEX_SRC})
target_include_directories(libqcex PUBLIC
    include
    ../common/include
)
target_link_libraries(libqcex
    PUBLIC
    libmatrix
    blasbind
    cint
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

# Apply build-specific flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(libqcex PRIVATE -Og -g3 -fno-openmp)
else()
    target_compile_options(libqcex PRIVATE -O3 -march=native)
    if(OpenMP_CXX_FOUND)
        target_compile_options(libqcex PRIVATE ${OpenMP_CXX_FLAGS})
        target_link_libraries(libqcex PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

# === Build each test in tests/*.cpp ===
file(GLOB LIBQCEX_TESTS tests/*.cpp)

foreach(test_file ${LIBQCEX_TESTS})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})

    target_include_directories(${test_name} PRIVATE
        include
        ../common/include
    )

    target_link_libraries(${test_name}
        PRIVATE
        libqcex
        libmatrix
        blasbind
        cint
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${test_name} PRIVATE -Og -g3 -fno-openmp)
    else()
        target_compile_options(${test_name} PRIVATE -O3 -march=native)
        if(OpenMP_CXX_FOUND)
            target_compile_options(${test_name} PRIVATE ${OpenMP_CXX_FLAGS})
            target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endif()
endforeach()