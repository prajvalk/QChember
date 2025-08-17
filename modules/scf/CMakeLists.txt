project(newscf_scf)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(newscf_scf SHARED scf.cpp)
target_link_libraries(newscf_scf newscf_libcint newscf_utils lacislapack lacis ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(he_hfm_test tests/he_hfm.cpp)
target_link_libraries(he_hfm_test newscf_scf newscf_utils newscf_libcint cint ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} lacis lacislapack)
add_dependencies(he_hfm_test copy_test_data)
target_compile_definitions(he_hfm_test PRIVATE
        NEWSCF_TEST_DIR="${TEST_DATA_DEST_DIR}"
)

enable_testing()
add_test(NAME helium.hfm.matrix.test COMMAND he_hfm_test)