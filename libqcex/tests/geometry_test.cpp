#include <fstream>
#include "testing_api.hpp"
#include "qcex.hpp"

using namespace newscf::testing;
using namespace newscf::qcex;

int main() {

    TestHandle hnd;
    hnd.testsuite = "qcex_geometry";

    init_tests();

    std::string water_xyz = "3 \n Water \n H    0.749023100410    0.000000000000    0.539419524871 \n O    0.000000000000    0.000000000000   -0.067976667956 \n H   -0.749023100410   -0.000000000000    0.539419524871";

    std::ofstream water_file ("test_water.xyz");
    water_file << water_xyz;
    water_file.close();

    add_test("water_test", hnd);

    double* WORK;

    load_molecule("test_water.xyz", &WORK);

    TEST_ASSERT_EQ (WORK[0], 1.0, hnd);
    TEST_ASSERT_EQ (WORK[1], 0.749023100410, hnd);
    TEST_ASSERT_EQ (WORK[2], 0.0, hnd);
    TEST_ASSERT_EQ (WORK[3], 0.539419524871, hnd);
    TEST_ASSERT_EQ (WORK[4], 8.0, hnd);
    TEST_ASSERT_EQ (WORK[5], 0.0, hnd);
    TEST_ASSERT_EQ (WORK[6], 0.0, hnd);
    TEST_ASSERT_EQ (WORK[7], -0.067976667956, hnd);
    TEST_ASSERT_EQ (WORK[8], 1.0, hnd);
    TEST_ASSERT_EQ (WORK[9], -0.749023100410, hnd);
    TEST_ASSERT_EQ (WORK[10], 0.0, hnd)
    TEST_ASSERT_EQ (WORK[11], 0.539419524871, hnd);

    destroy_molecule(WORK);

    complete_test (hnd);
    end_tests(hnd);

    return hnd.exitcode;

}