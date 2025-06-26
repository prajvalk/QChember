#include "libcint_wrapper.hpp"
#include "qcex.hpp"
#include "testing_api.hpp"

using namespace newscf::testing;
using namespace newscf::qcex;

int main() {
    TestHandle h;
    h.testsuite = "qcex_libcint_wrapper";

    init_tests();

    add_test("water/def2-TZVPP", h);

    double* geom;
    double* basis;

    load_molecule("test_water.xyz", &geom);
    load_basis("def2-TZVPP.bas", &basis);

    int* atm;
    int* bas;
    int natm, nbas;
    double* env;
    int nenv;

    prepare_libcint_data (geom, basis, &atm, &natm, &bas, &nbas, &env, &nenv);

    destroy_environment(atm, bas, env);
    destroy_basis(basis);
    destroy_molecule(geom);

    end_tests(h);

    return h.exitcode;
}