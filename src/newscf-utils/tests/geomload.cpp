#include "newscf/newscf.hpp"
#include "newscf/testing_cis.hpp"

int main () {
    using namespace newscf::cis;
    using namespace newscf;
    using namespace newscf::testing_cis;

    TestHandle hnd;
    hnd.testsuite = "qcex_geometry";

    init_tests();

    add_test("water_test", hnd);

    double* WORK;

    std::string str = NEWSCF_TEST_DIR;
    str += "/geometries/testsuite_bronze/01_Water.xyz";
    int mol_size = load_molecule(str, &WORK);

    TEST_ASSERT_EQ (mol_size, 3, hnd); // Number of Atoms
    TEST_ASSERT_EQ (WORK[0], 1.0, hnd); // Z-Val of H
    TEST_ASSERT_EQ (WORK[1], 0.752700863259, hnd); // X
    TEST_ASSERT_EQ (WORK[2], 0.0, hnd);            // Y
    TEST_ASSERT_EQ (WORK[3], 0.525364460618, hnd); // Z
    TEST_ASSERT_EQ (WORK[4], 8.0, hnd); // Z-Val of O
    TEST_ASSERT_EQ (WORK[5], 0.0, hnd); // X
    TEST_ASSERT_EQ (WORK[6], 0.0, hnd); // Y
    TEST_ASSERT_EQ (WORK[7], -0.066205474308, hnd); //Z
    TEST_ASSERT_EQ (WORK[8], 1.0, hnd); // Z-Val of H
    TEST_ASSERT_EQ (WORK[9], -0.752700863259, hnd); // X
    TEST_ASSERT_EQ (WORK[10], 0.0, hnd); // Y
    TEST_ASSERT_EQ (WORK[11], 0.525364460618, hnd); // Z

    destroy_molecule(WORK);

    complete_test (hnd);
    end_tests(hnd);

    return hnd.exitcode;
}