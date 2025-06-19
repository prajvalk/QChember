#ifndef NEWSCF_BLASBIND_HPP
#define NEWSCF_BLASBIND_HPP

#include "commons.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

/* SIMD Fallback Routines */
void blasbind_sgemm (const char TRANSA, const char TRANSB, const int M, const int N, const int K,
                const float ALPHA, const float* A, const int LDA,
                const float* B, const int LDB, const float BETA,
                float* C, const int LDC);
void blasbind_dgemm (const char TRANSA, const char TRANSB, const int M, const int N, const int K,
                const double ALPHA, const double* A, const int LDA,
                const double* B, const int LDB, const double BETA,
                double* C, const int LDC);
template <typename T>
void blasbind_xgemm (const char TRANSA, const char TRANSB, const int M, const int N, const int K,
                const double ALPHA, const double* A, const int LDA,
                const double* B, const int LDB, const double BETA,
                double* C, const int LDC) {
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            T total = T(0);
            for (int k = 0; k < K; ++k) {
                T a_val = (TRANSA == 'T' || TRANSA == 't') ? A[k + i * LDA] : A[i + k * LDA];
                T b_val = (TRANSB == 'T' || TRANSB == 't') ? B[j + k * LDB] : B[k + j * LDB];
                total += a_val * b_val;
            }
            C[i + j * LDC] = ALPHA * total + BETA * C[i + j * LDC];
        }
    }
}

namespace newscf::blasbind {

    template <typename T>
    inline void matmul(const char TRANSA, const char TRANSB, 
                   const size_t M, const size_t N, const size_t K,
                   const T ALPHA, const T* A, size_t LDA, 
                   const T* B, size_t LDB, const T BETA, 
                   T* C, size_t LDC) {
        // Compute default leading dimensions if -1 was passed
        const bool transA = (TRANSA == 'T' || TRANSA == 't');
        const bool transB = (TRANSB == 'T' || TRANSB == 't');

        size_t lda = (LDA == size_t(-1)) ? (transA ? K : M) : LDA;
        size_t ldb = (LDB == size_t(-1)) ? (transB ? N : K) : LDB;
        size_t ldc = (LDC == size_t(-1)) ? M : LDC;

        if constexpr (std::is_same<T, float>::value) {
            blasbind_sgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, lda, B, ldb, BETA, C, ldc);
        } else if constexpr (std::is_same<T, double>::value) {
            blasbind_dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, lda, B, ldb, BETA, C, ldc);
        } else {
            blasbind_xgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, lda, B, ldb, BETA, C, ldc);
        }
    }

    template <typename T>
    inline void matmul(const size_t M, const size_t N, const size_t K,
                    const T* A, const T* B, T* C,
                    const char TRANSA = 'N', const char TRANSB = 'N',
                    const T ALPHA = T(1), const T BETA = T(0),
                    const size_t LDA = size_t(-1), const size_t LDB = size_t(-1), const size_t LDC = size_t(-1)) {
        matmul<T>(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
    }


}

#endif