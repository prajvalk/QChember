project(blasbind LANGUAGES CXX)

file(GLOB_RECURSE BLASBIND_SRC src/*.cpp)

add_library(blasbind SHARED ${BLASBIND_SRC})
target_include_directories(blasbind PUBLIC include ../common/include)
target_link_libraries(blasbind PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Set target-specific flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(blasbind PRIVATE -Og -g3 -fno-openmp)
else()
    target_compile_options(blasbind PRIVATE -O3 -march=native)
    if(OpenMP_CXX_FOUND)
        target_compile_options(blasbind PRIVATE ${OpenMP_CXX_FLAGS})
        target_link_libraries(blasbind PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

# Build all tests
file(GLOB BLASBIND_TESTS tests/*.cpp)
foreach(test_file ${BLASBIND_TESTS})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} PRIVATE blasbind)
    target_include_directories(${test_name} PRIVATE include ../common/include)

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
endforeach()
