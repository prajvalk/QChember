#include "blasbind.hpp"

#include <immintrin.h>
#include <omp.h>

/* Export Bindings to FBLAS */
#if defined(NEWSCF_BLAS_LINK)
extern "C" {
    void sgemm_ (char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
}
#endif

void blasbind_sgemm (const char TRANSA, const char TRANSB, const int M, const int N, const int K,
                const float ALPHA, const float* A, const int LDA,
                const float* B, const int LDB, const float BETA,
                float* C, const int LDC) {
#if defined(NEWSCF_BLAS_LINK)
    sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
    return;
#endif

    const bool transA = (TRANSA == 'T' || TRANSA == 't');
    const bool transB = (TRANSB == 'T' || TRANSB == 't');

#if defined(__AVX512F__)
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            __m512 sum = _mm512_setzero_ps();
            int k = 0;
            for (; k + 15 < K; k += 16) {
                float a_buf[16], b_buf[16];
                for (int z = 0; z < 16; ++z) {
                    a_buf[z] = transA ? A[k + z + i * LDA] : A[i + (k + z) * LDA]; // A(i,k)
                    b_buf[z] = transB ? B[j + (k + z) * LDB] : B[(k + z) + j * LDB]; // B(k,j)
                }
                __m512 a = _mm512_loadu_ps(a_buf);
                __m512 b = _mm512_loadu_ps(b_buf);
                sum = _mm512_fmadd_ps(a, b, sum);
            }
            float buffer[16];
            _mm512_storeu_ps(buffer, sum);
            float total = 0.0f;
            for (int z = 0; z < 16; ++z) total += buffer[z];
            for (; k < K; ++k) {
                float a_val = transA ? A[k + i * LDA] : A[i + k * LDA];
                float b_val = transB ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }

#elif defined(__AVX2__)
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            __m256 sum = _mm256_setzero_ps();
            int k = 0;
            for (; k + 7 < K; k += 8) {
                float a_buf[8], b_buf[8];
                for (int z = 0; z < 8; ++z) {
                    a_buf[z] = transA ? A[k + z + i * LDA] : A[i + (k + z) * LDA];
                    b_buf[z] = transB ? B[j + (k + z) * LDB] : B[(k + z) + j * LDB];
                }
                __m256 a = _mm256_loadu_ps(a_buf);
                __m256 b = _mm256_loadu_ps(b_buf);
                sum = _mm256_fmadd_ps(a, b, sum);
            }
            float buffer[8];
            _mm256_storeu_ps(buffer, sum);
            float total = 0.0f;
            for (int z = 0; z < 8; ++z) total += buffer[z];
            for (; k < K; ++k) {
                float a_val = transA ? A[k + i * LDA] : A[i + k * LDA];
                float b_val = transB ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }

#else
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            float total = 0.0f;
            for (int k = 0; k < K; ++k) {
                float a_val = transA ? A[k + i * LDA] : A[i + k * LDA];
                float b_val = transB ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }
#endif
}