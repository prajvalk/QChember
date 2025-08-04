#include "newscf/cis.hpp"
#include "newscf/testing_cis.hpp"

template<newscf::cis::GeneralizedEigenvalueSolverType type>
void run_geigensolver_test (const std::string& name, newscf::testing_cis::TestHandle& handle) {
    using namespace newscf::cis;
    using namespace newscf::testing_cis;

    add_test(name + "Problem #1", handle);

    LACIS<DENSE_LAPACK, double> matA;
    matA.resizeToMatrix(3, 3);

    LACIS<DENSE_LAPACK, double> matB;
    matB.resizeToMatrix(3, 3);

    matA.matrixSet(0, 0, 4);
    matA.matrixSet(0, 1, 1);
    matA.matrixSet(0, 2, 1);
    matA.matrixSet(1, 0, 1);
    matA.matrixSet(2, 0, 1);
    matA.matrixSet(1, 1, 3);
    matA.matrixSet(2, 2, 2);
    matA.matrixSet(1, 2, 0);
    matA.matrixSet(2, 1, 0);

    matB.matrixSet(0, 0, 3);
    matB.matrixSet(0, 1, 1);
    matB.matrixSet(0, 2, 0);
    matB.matrixSet(1, 0, 1);
    matB.matrixSet(2, 0, 0);
    matB.matrixSet(1, 1, 2);
    matB.matrixSet(2, 2, 1);
    matB.matrixSet(1, 2, 0);
    matB.matrixSet(2, 1, 0);

    GeneralizedEigenvalueSolver<type> solver;
    solver.set_itype(1);
    solver.set_jobz('V');
    solver.set_uplo('U');
    solver.set_n(3);
    solver.init();

    solver.set_a(matA);
    solver.set_b(matB);

    solver.solve();

    LACIS<DENSE_LAPACK, double> sol_eval;
    sol_eval.resizeToVector(3);
    LACIS<DENSE_LAPACK, double> sol_evec;
    sol_evec.resizeToMatrix(3, 3);

    solver.copy_eigenvalues(sol_eval);
    solver.copy_eigenvectors(sol_evec);

    TEST_ASSERT_EQ_TOL(sol_eval.vectorGet(0), static_cast<double>(1), 1e-6, handle);
    TEST_ASSERT_EQ_TOL(sol_eval.vectorGet(1), static_cast<double>(1.5527864), 1e-6, handle);
    TEST_ASSERT_EQ_TOL(sol_eval.vectorGet(2), static_cast<double>(2.4472136), 1e-6, handle);

    complete_test(handle);
}

int main() {
    using namespace newscf::cis;
    using namespace newscf::testing_cis;

    TestHandle hnd;
    hnd.testsuite = "lacis-lapack tests";

    init_tests();

    run_geigensolver_test<DSYGV>("DSYGV", hnd);

    end_tests(hnd);

    return hnd.exitcode;
}