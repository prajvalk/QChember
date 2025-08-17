project(newscf_utils)

add_library(newscf_utils SHARED geometry.cpp bse_handler.cpp)

add_executable(geomload_test tests/geomload.cpp)
target_link_libraries(geomload_test newscf_utils)
add_dependencies(geomload_test copy_test_data)
target_compile_definitions(geomload_test PRIVATE
        NEWSCF_TEST_DIR="${TEST_DATA_DEST_DIR}"
)

add_executable(bseload_test tests/bseload.cpp)
target_link_libraries(bseload_test newscf_utils)
add_dependencies(bseload_test copy_test_data)
target_compile_definitions(bseload_test PRIVATE
        NEWSCF_TEST_DIR="${TEST_DATA_DEST_DIR}"
)

enable_testing()
add_test(NAME geomload.test COMMAND geomload_test)
add_test(NAME bseload.test COMMAND bseload_test)