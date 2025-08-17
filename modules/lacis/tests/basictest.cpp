#include "newscf/cis.hpp"
#include "newscf/testing_cis.hpp"

int main() {

    using namespace newscf::cis;
    using namespace newscf::testing_cis;

    TestHandle test;
    test.testsuite = "lacis basic testsuite";
    init_tests();

    add_test("Vector Operations", test);

    LACIS<DENSE_LAPACK, double> vector;
    vector.resizeToVector(5);

    for (int i = 0; i < 5; i++) vector.vectorSet(i, 2 * i);

    for (int i = 0; i < 5; i++)
        TEST_ASSERT_EQ(vector.vectorGet(i), static_cast<double>(2 * i), test);

    complete_test(test);

    add_test("Matrix Operations", test);

    LACIS<DENSE_LAPACK, double> matrix;
    matrix.resizeToMatrix(10, 10);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            matrix.matrixSet(i, j, i * j);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            TEST_ASSERT_EQ(matrix.matrixGet(i, j), static_cast<double>(i * j), test);

    complete_test(test);

    add_test("Tensor3D Operations", test);

    LACIS<DENSE_LAPACK, double> tensor3D;
    tensor3D.resizeToTensor3D(10, 10, 10);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                tensor3D.tensor3DSet(i, j, k, i * j * k);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                TEST_ASSERT_EQ(tensor3D.tensor3DGet(i, j, k), static_cast<double>(i * j * k), test);

    complete_test(test);

    add_test("Tensor4D Operations", test);

    LACIS<DENSE_LAPACK, double> tensor4D;
    tensor4D.resizeToTensor4D(10, 10, 10, 10);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                for (int l = 0; l < 10; l++)
                    tensor4D.tensor4DSet(i, j, k, l, i * j + k * l);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                for (int l = 0; l < 10; l++)
                    TEST_ASSERT_EQ(tensor4D.tensor4DGet(i, j, k, l), static_cast<double>(i * j + k * l), test);

    complete_test(test);

    add_test("Copy Test", test);

    LACIS<DENSE_LAPACK, double> copy;
    copy = tensor4D;
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                for (int l = 0; l < 10; l++)
                    TEST_ASSERT_EQ(copy.tensor4DGet(i, j, k, l), static_cast<double>(i * j + k * l), test);

    complete_test(test);

    end_tests(test);
    return test.exitcode;
}