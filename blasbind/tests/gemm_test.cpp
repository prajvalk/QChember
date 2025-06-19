#include "testing_api.hpp"
#include "blasbind.hpp"

using namespace newscf::testing;

int main() {
    TestHandle h;
    h.testsuite = "newgemm::blasbind::gemm";
    init_tests();

    constexpr int M = 2, N = 2, K = 2;
    const float A[M * K] = {
        1.0f, 2.0f,
        3.0f, 4.0f
    }; // Column-major

    const float B[K * N] = {
        5.0f, 6.0f,
        7.0f, 8.0f
    };

    float C[M * N] = {
        0.0f, 0.0f,
        0.0f, 0.0f
    };

    const float expected[M * N] = {
        23.0f, 34.0f,
        31.0f, 46.0f
    };

    const float alpha = 1.0f;
    const float beta = 0.0f;

    add_test ("blasbind_sgemm", h);

    blasbind_sgemm('N', 'N', M, N, K, alpha, A, M, B, K, beta, C, M);

    for (size_t i = 0; i < M * N; i++)
        TEST_ASSERT_EQ (C[i], expected[i], h);

    end_tests(h);
}