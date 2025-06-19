#include "blasbind.hpp"

#include <immintrin.h>
#include <omp.h>

/* Export Bindings to FBLAS */
#if defined(NEWSCF_BLAS_LINK)
extern "C" {
    void dgemm_ (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
}
#endif

void blasbind_dgemm(const char TRANSA, const char TRANSB, const int M, const int N, const int K,
                    const double ALPHA, const double* A, const int LDA,
                    const double* B, const int LDB, const double BETA,
                    double* C, const int LDC) {
#if defined(NEWSCF_BLAS_LINK)
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
    return;
#endif

    const bool transA = (TRANSA == 'T' || TRANSA == 't');
    const bool transB = (TRANSB == 'T' || TRANSB == 't');

#if defined(__AVX512F__)
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            __m512d sum = _mm512_setzero_pd();
            int k = 0;
            for (; k + 7 < K; k += 8) {
                double a_buf[8], b_buf[8];
                for (int z = 0; z < 8; ++z) {
                    a_buf[z] = transA ? A[k + z + i * LDA] : A[i + (k + z) * LDA];  // A(i,k)
                    b_buf[z] = transB ? B[j + (k + z) * LDB] : B[(k + z) + j * LDB];  // B(k,j)
                }
                __m512d a = _mm512_loadu_pd(a_buf);
                __m512d b = _mm512_loadu_pd(b_buf);
                sum = _mm512_fmadd_pd(a, b, sum);
            }
            double buffer[8];
            _mm512_storeu_pd(buffer, sum);
            double total = 0.0;
            for (int z = 0; z < 8; ++z) total += buffer[z];
            for (; k < K; ++k) {
                double a_val = transA ? A[k + i * LDA] : A[i + k * LDA];
                double b_val = transB ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }

#elif defined(__AVX2__)
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            __m256d sum = _mm256_setzero_pd();
            int k = 0;
            for (; k + 3 < K; k += 4) {
                double a_buf[4], b_buf[4];
                for (int z = 0; z < 4; ++z) {
                    a_buf[z] = transA ? A[k + z + i * LDA] : A[i + (k + z) * LDA];
                    b_buf[z] = transB ? B[j + (k + z) * LDB] : B[(k + z) + j * LDB];
                }
                __m256d a = _mm256_loadu_pd(a_buf);
                __m256d b = _mm256_loadu_pd(b_buf);
                sum = _mm256_fmadd_pd(a, b, sum);
            }
            double buffer[4];
            _mm256_storeu_pd(buffer, sum);
            double total = 0.0;
            for (int z = 0; z < 4; ++z) total += buffer[z];
            for (; k < K; ++k) {
                double a_val = transA ? A[k + i * LDA] : A[i + k * LDA];
                double b_val = transB ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }

#else
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            double total = 0.0;
            for (int k = 0; k < K; ++k) {
                double a_val = transA ? A[k + i * LDA] : A[i + k * LDA];
                double b_val = transB ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }
#endif
}
