#include "matrix.hpp"

#include <immintrin.h>
#include <omp.h>
#include <type_traits>
#include <stdexcept>

namespace newscf::libmatrix {

    // ======== FLOAT ADDITION: C = A + B ========
    void add_float_simd(const float* A, const float* B, float* C, size_t N) {
    #if defined(__AVX512F__)
        constexpr size_t simd_width = 16;
    #elif defined(__AVX2__)
        constexpr size_t simd_width = 8;
    #else
        constexpr size_t simd_width = 1;
    #endif

        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            if constexpr (simd_width > 1) {
                if (i + simd_width <= N) {
                #if defined(__AVX512F__)
                    __m512 a = _mm512_loadu_ps(A + i);
                    __m512 b = _mm512_loadu_ps(B + i);
                    __m512 c = _mm512_add_ps(a, b);
                    _mm512_storeu_ps(C + i, c);
                #elif defined(__AVX2__)
                    __m256 a = _mm256_loadu_ps(A + i);
                    __m256 b = _mm256_loadu_ps(B + i);
                    __m256 c = _mm256_add_ps(a, b);
                    _mm256_storeu_ps(C + i, c);
                #endif
                    i += simd_width - 1;  // skip next simd_width-1 elements
                } else {
                    C[i] = A[i] + B[i];
                }
            } else {
                C[i] = A[i] + B[i];
            }
        }
    }

    // ======== DOUBLE ADDITION: C = A + B ========
    void add_double_simd(const double* A, const double* B, double* C, size_t N) {
    #if defined(__AVX512F__)
        constexpr size_t simd_width = 8;
    #elif defined(__AVX2__)
        constexpr size_t simd_width = 4;
    #else
        constexpr size_t simd_width = 1;
    #endif

        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            if constexpr (simd_width > 1) {
                if (i + simd_width <= N) {
                #if defined(__AVX512F__)
                    __m512d a = _mm512_loadu_pd(A + i);
                    __m512d b = _mm512_loadu_pd(B + i);
                    __m512d c = _mm512_add_pd(a, b);
                    _mm512_storeu_pd(C + i, c);
                #elif defined(__AVX2__)
                    __m256d a = _mm256_loadu_pd(A + i);
                    __m256d b = _mm256_loadu_pd(B + i);
                    __m256d c = _mm256_add_pd(a, b);
                    _mm256_storeu_pd(C + i, c);
                #endif
                    i += simd_width - 1;
                } else {
                    C[i] = A[i] + B[i];
                }
            } else {
                C[i] = A[i] + B[i];
            }
        }
    }

    // ======== IN-PLACE FLOAT ADDITION: C += A ========
    void add_inplace_float_simd(float* C, const float* A, size_t N) {
    #if defined(__AVX512F__)
        constexpr size_t simd_width = 16;
    #elif defined(__AVX2__)
        constexpr size_t simd_width = 8;
    #else
        constexpr size_t simd_width = 1;
    #endif

        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            if constexpr (simd_width > 1) {
                if (i + simd_width <= N) {
                #if defined(__AVX512F__)
                    __m512 a = _mm512_loadu_ps(A + i);
                    __m512 c = _mm512_loadu_ps(C + i);
                    c = _mm512_add_ps(a, c);
                    _mm512_storeu_ps(C + i, c);
                #elif defined(__AVX2__)
                    __m256 a = _mm256_loadu_ps(A + i);
                    __m256 c = _mm256_loadu_ps(C + i);
                    c = _mm256_add_ps(a, c);
                    _mm256_storeu_ps(C + i, c);
                #endif
                    i += simd_width - 1;
                } else {
                    C[i] += A[i];
                }
            } else {
                C[i] += A[i];
            }
        }
    }

    // ======== IN-PLACE DOUBLE ADDITION: C += A ========
    void add_inplace_double_simd(double* C, const double* A, size_t N) {
    #if defined(__AVX512F__)
        constexpr size_t simd_width = 8;
    #elif defined(__AVX2__)
        constexpr size_t simd_width = 4;
    #else
        constexpr size_t simd_width = 1;
    #endif

        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            if constexpr (simd_width > 1) {
                if (i + simd_width <= N) {
                #if defined(__AVX512F__)
                    __m512d a = _mm512_loadu_pd(A + i);
                    __m512d c = _mm512_loadu_pd(C + i);
                    c = _mm512_add_pd(a, c);
                    _mm512_storeu_pd(C + i, c);
                #elif defined(__AVX2__)
                    __m256d a = _mm256_loadu_pd(A + i);
                    __m256d c = _mm256_loadu_pd(C + i);
                    c = _mm256_add_pd(a, c);
                    _mm256_storeu_pd(C + i, c);
                #endif
                    i += simd_width - 1;
                } else {
                    C[i] += A[i];
                }
            } else {
                C[i] += A[i];
            }
        }
    }

}  // namespace newscf::libmatrix
