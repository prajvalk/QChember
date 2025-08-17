#include "newscf/newscf.hpp"
#include "newscf/testing_cis.hpp"

int main() {
    using namespace newscf;
    using namespace newscf::cis;
    using namespace newscf::testing_cis;

    TestHandle h;
    h.testsuite = "Helium HF/def2-SVP test";
    init_tests();

    double* mol;
    double* bas;

    std::string testdir = NEWSCF_TEST_DIR;

    const int molsize = load_molecule(testdir + "/geometries/testsuite_bronze/11_Helium.xyz", &mol);
    const int bassize = load_basis(testdir + "/basis_sets/def2-svp.1.gbs", &bas, mol, molsize);

    IntegralEngineHandle<LIBCINT> ih;
    initialize_handle(mol, molsize, bas, bassize, ih);

    double en = run_rhf<LIBCINT, DENSE_LAPACK>(ih);

    TEST_ASSERT_EQ_TOL(en, -2.855160, 1e-6, h);

    destroy_handle(ih);
    destroy_molecule(mol);
    destroy_basis(bas);

    complete_test(h);

    end_tests(h);
    return h.exitcode;
}