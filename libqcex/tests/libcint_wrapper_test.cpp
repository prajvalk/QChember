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

    double* geom;
    double* basis;

    int molecule_size = load_molecule("test_water.xyz", &geom);
    int basis_size = load_basis("def2-SVP.bas", &basis);

    IntegralEngineHandle handle;

    intialize_handle(geom, molecule_size, basis, basis_size, &handle);
    Matrix<double>* ovp = nullptr;
    calculate_overlap_matrix (&handle, &ovp);

    for (int i = 0; i < ovp->sz_rows; i++) {
        for (int j = 0; j < ovp->sz_cols; j++) {
            if(ovp->get(i, j) != 0) std::cout << i << ", " << j << ", " << ovp->get(i, j) << "\n";
        }
    }

    destroy_handle(&handle);
    destroy_basis(basis);
    destroy_molecule(geom);

    end_tests(h);

    return h.exitcode;
}