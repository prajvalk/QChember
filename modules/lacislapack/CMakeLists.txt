project(lacislapack)

add_library(lacislapack SHARED symlinearsolver.cpp symeigen.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(lacislapack_linearsolver_test tests/linearsolver_test.cpp)
target_link_libraries(lacislapack_linearsolver_test lacislapack lacis ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(lacislapack_eigensolver_test tests/eigensolver_test.cpp)
target_link_libraries(lacislapack_eigensolver_test lacislapack lacis ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

enable_testing()
add_test(NAME linearsolver.test COMMAND lacislapack_linearsolver_test)
add_test(NAME eigensolver.test COMMAND lacislapack_eigensolver_test)