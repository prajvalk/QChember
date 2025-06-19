#include "matrix.hpp"
#include "testing_api.hpp"

using namespace newscf::libmatrix;
using namespace newscf::testing;

template<typename T>
void run_matrix_tests(TestHandle& h, const std::string& prefix) {
    size_t rows = 1024;
    size_t cols = 2048;

    Matrix<T> A(rows, cols);
    Matrix<T> B(rows, cols);
    Matrix<T> C, D;

    // Initialize A(i,j) = i + j, B(i,j) = j - i
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++) {
            A.set(i, j, T(i) + T(j));
            B.set(i, j, T(j) - T(i));
        }

    // Addition
    add_test(prefix + "::add", h);
    C = A + B;
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++)
            TEST_ASSERT_EQ (C.get(i, j), T(2 * j), h);

    // In-place addition
    add_test(prefix + "::add_inplace", h);
    D = A;
    D += B;
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++)
            TEST_ASSERT_EQ (D.get(i, j), T(2 * j), h);

    // Subtraction
    add_test(prefix + "::sub", h);
    C = A - B;
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++)
            TEST_ASSERT_EQ (C.get(i, j), T(2 * i), h);

    // In-place subtraction
    add_test(prefix + "::sub_inplace", h);
    D = A;
    D -= B;
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++)
            TEST_ASSERT_EQ (D.get(i, j), T(2 * i), h);

// Scalar multiplication
    T alpha = T(3);
    add_test(prefix + "::scalar_mul", h);
    C = A * alpha;
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++)
            TEST_ASSERT_EQ (C.get(i, j), T(3 * (i + j)), h);

    // In-place scalar multiplication
    add_test(prefix + "::scalar_mul_inplace", h);
    D = A;
    D *= alpha;
    for (size_t j = 0; j < cols; j++)
        for (size_t i = 0; i < rows; i++)
            TEST_ASSERT_EQ (D.get(i, j), T(3 * (i + j)), h);
}


int main() {
    TestHandle h;
    h.testsuite = "newgemm::libmatrix";
    init_tests();

    run_matrix_tests<float>(h, "float");
    run_matrix_tests<double>(h, "double");
    run_matrix_tests<long>(h, "long");

    end_tests(h);
    return h.exitcode;
}