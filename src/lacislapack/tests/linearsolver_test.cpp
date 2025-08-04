#include "newscf/cis.hpp"
#include "newscf/testing_cis.hpp"

template <newscf::cis::LinearSolverType T>
void run_linear_solver_test(const std::string& str, newscf::testing_cis::TestHandle& handle) {
    using namespace newscf::testing_cis;
    using namespace newscf::cis;

    add_test(str+" Problem #1", handle);

    LACIS<DENSE_LAPACK, double> matA;
    matA.resizeToMatrix(3, 3);

    matA.matrixSet(0, 0, 4);
    matA.matrixSet(0, 1, 1);
    matA.matrixSet(0, 2, 2);
    matA.matrixSet(1, 0, 1);
    matA.matrixSet(1, 1, 3);
    matA.matrixSet(1, 2, 0);
    matA.matrixSet(2, 0, 2);
    matA.matrixSet(2, 1, 0);
    matA.matrixSet(2, 2, 5);

    LACIS<DENSE_LAPACK, double> vecB;
    vecB.resizeToVector(3);

    vecB.vectorSet(0, 7);
    vecB.vectorSet(1, 8);
    vecB.vectorSet(2, 9);

    LACIS<DENSE_LAPACK, double> vecX;
    vecX.resizeToVector(3);

    LinearSolver<T> solver;

    solver.set_n(3);
    solver.set_nrhs(1);
    solver.set_uplo('U');
    solver.init();

    solver.set_a(matA);
    solver.set_b(vecB);

    solver.solve();

    solver.copy_result(vecX);

    TEST_ASSERT_EQ_TOL(vecX.vectorGet(0), 0.25581395, 1e-6, handle);
    TEST_ASSERT_EQ_TOL(vecX.vectorGet(1), 2.58139535, 1e-6, handle);
    TEST_ASSERT_EQ_TOL(vecX.vectorGet(2), 1.69767442, 1e-6, handle);

    complete_test(handle);
}

int main() {
    using namespace newscf::cis;
    using namespace newscf::testing_cis;

    TestHandle hnd;
    hnd.testsuite = "lacis-lapack tests";

    init_tests();

    run_linear_solver_test<DSYSV>("DSYSV", hnd);
    run_linear_solver_test<DSYSV_ROOK>("DSYSV_ROOK", hnd);
    run_linear_solver_test<DSYSV_RK>("DSYSV_RK", hnd);
    run_linear_solver_test<DSYSV_AA>("DSYSV_AA", hnd);
    run_linear_solver_test<DSYSV_AA_2>("DSYSV_AA_2", hnd);

    end_tests(hnd);

    return hnd.exitcode;
}