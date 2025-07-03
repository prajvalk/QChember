#include "qcex.hpp"
#include "testing_api.hpp"
#include <cstring>

using namespace newscf::testing;
using namespace newscf::qcex;

int main() {
    TestHandle h;
    h.testsuite = "qcex_libcint_wrapper";

    init_tests();

    add_test("water/def2-TZVPP", h);

    double* geom = new double[4] {2, 0, 0, 0};
    double* basis;

    int molecule_size = 1;
    int basis_size = load_basis("def2-SVP.bas", &basis);

    IntegralEngineHandle handle;

    intialize_handle(geom, molecule_size, basis, basis_size, &handle);
    
    double* ener = nullptr;
    atom_rhf_energies (&handle, &ener);

    destroy_handle(&handle);
    destroy_basis(basis);
    destroy_molecule(geom);

    end_tests(h);

    return h.exitcode;
}