#include "matrix.hpp"

#include <immintrin.h>
#include <omp.h>

namespace newscf::libmatrix {

    // -------- FLOAT SCALAR MULTIPLICATION --------
    void mul_scalar_float_simd(const float* A, float* C, float scalar, size_t N) {
    #if defined(__AVX512F__)
        __m512 s = _mm512_set1_ps(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 16) {
            if (i + 16 <= N) {
                __m512 a = _mm512_loadu_ps(A + i);
                __m512 c = _mm512_mul_ps(a, s);
                _mm512_storeu_ps(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] = A[j] * scalar;
            }
        }

    #elif defined(__AVX2__)
        __m256 s = _mm256_set1_ps(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 8) {
            if (i + 8 <= N) {
                __m256 a = _mm256_loadu_ps(A + i);
                __m256 c = _mm256_mul_ps(a, s);
                _mm256_storeu_ps(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] = A[j] * scalar;
            }
        }

    #else
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] = A[i] * scalar;
    #endif
    }

    // -------- DOUBLE SCALAR MULTIPLICATION --------
    void mul_scalar_double_simd(const double* A, double* C, double scalar, size_t N) {
    #if defined(__AVX512F__)
        __m512d s = _mm512_set1_pd(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 8) {
            if (i + 8 <= N) {
                __m512d a = _mm512_loadu_pd(A + i);
                __m512d c = _mm512_mul_pd(a, s);
                _mm512_storeu_pd(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] = A[j] * scalar;
            }
        }

    #elif defined(__AVX2__)
        __m256d s = _mm256_set1_pd(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 4) {
            if (i + 4 <= N) {
                __m256d a = _mm256_loadu_pd(A + i);
                __m256d c = _mm256_mul_pd(a, s);
                _mm256_storeu_pd(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] = A[j] * scalar;
            }
        }

    #else
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] = A[i] * scalar;
    #endif
    }

    // -------- INPLACE FLOAT SCALAR MULTIPLICATION --------
    void mul_scalar_inplace_float_simd(float* C, float scalar, size_t N) {
    #if defined(__AVX512F__)
        __m512 s = _mm512_set1_ps(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 16) {
            if (i + 16 <= N) {
                __m512 c = _mm512_loadu_ps(C + i);
                c = _mm512_mul_ps(c, s);
                _mm512_storeu_ps(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] *= scalar;
            }
        }

    #elif defined(__AVX2__)
        __m256 s = _mm256_set1_ps(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 8) {
            if (i + 8 <= N) {
                __m256 c = _mm256_loadu_ps(C + i);
                c = _mm256_mul_ps(c, s);
                _mm256_storeu_ps(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] *= scalar;
            }
        }

    #else
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] *= scalar;
    #endif
    }

    // -------- INPLACE DOUBLE SCALAR MULTIPLICATION --------
    void mul_scalar_inplace_double_simd(double* C, double scalar, size_t N) {
    #if defined(__AVX512F__)
        __m512d s = _mm512_set1_pd(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 8) {
            if (i + 8 <= N) {
                __m512d c = _mm512_loadu_pd(C + i);
                c = _mm512_mul_pd(c, s);
                _mm512_storeu_pd(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] *= scalar;
            }
        }

    #elif defined(__AVX2__)
        __m256d s = _mm256_set1_pd(scalar);
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 4) {
            if (i + 4 <= N) {
                __m256d c = _mm256_loadu_pd(C + i);
                c = _mm256_mul_pd(c, s);
                _mm256_storeu_pd(C + i, c);
            } else {
                for (size_t j = i; j < N; ++j)
                    C[j] *= scalar;
            }
        }

    #else
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] *= scalar;
    #endif
    }

}

/*
// ------ MATRIX - MATRIX MULTIPLICATION ---------
extern "C" {
    void sgemm_(const char*, const char*, const int*, const int*, const int*,
                const float*, const float*, const int*, const float*, const int*,
                const float*, float*, const int*);

    void dgemm_(const char*, const char*, const int*, const int*, const int*,
                const double*, const double*, const int*, const double*, const int*,
                const double*, double*, const int*);
}

// Matrix C = A * B (column-major)
template <typename T>
void matmul_simd_fallback(const T* A, const T* B, T* C, size_t M, size_t N, size_t K);

// Float specialization
template <>
void matmul_simd_fallback<float>(const float* A, const float* B, float* C, size_t M, size_t N, size_t K) {
#if defined(__AVX512F__)
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            __m512 sum = _mm512_setzero_ps();
            size_t k = 0;
            for (; k + 15 < K; k += 16) {
                __m512 a = _mm512_loadu_ps(A + i + k * M);
                __m512 b = _mm512_loadu_ps(B + k + j * K);
                sum = _mm512_fmadd_ps(a, b, sum);
            }
            float buffer[16];
            _mm512_storeu_ps(buffer, sum);
            float total = 0.0f;
            for (int z = 0; z < 16; ++z) total += buffer[z];
            for (; k < K; ++k)
                total += A[i + k * M] * B[k + j * K];
            C[i + j * M] = total;
        }
    }
#elif defined(__AVX2__)
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            __m256 sum = _mm256_setzero_ps();
            size_t k = 0;
            for (; k + 7 < K; k += 8) {
                __m256 a = _mm256_loadu_ps(A + i + k * M);
                __m256 b = _mm256_loadu_ps(B + k + j * K);
                sum = _mm256_fmadd_ps(a, b, sum);
            }
            float buffer[8];
            _mm256_storeu_ps(buffer, sum);
            float total = 0.0f;
            for (int z = 0; z < 8; ++z) total += buffer[z];
            for (; k < K; ++k)
                total += A[i + k * M] * B[k + j * K];
            C[i + j * M] = total;
        }
    }
#else
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N; ++j)
        for (size_t i = 0; i < M; ++i) {
            float sum = 0;
            for (size_t k = 0; k < K; ++k)
                sum += A[i + k * M] * B[k + j * K];
            C[i + j * M] = sum;
        }
#endif
}

// Double specialization
template <>
void matmul_simd_fallback<double>(const double* A, const double* B, double* C, size_t M, size_t N, size_t K) {
#if defined(__AVX512F__)
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            __m512d sum = _mm512_setzero_pd();
            size_t k = 0;
            for (; k + 7 < K; k += 8) {
                __m512d a = _mm512_loadu_pd(A + i + k * M);
                __m512d b = _mm512_loadu_pd(B + k + j * K);
                sum = _mm512_fmadd_pd(a, b, sum);
            }
            double buffer[8];
            _mm512_storeu_pd(buffer, sum);
            double total = 0.0;
            for (int z = 0; z < 8; ++z) total += buffer[z];
            for (; k < K; ++k)
                total += A[i + k * M] * B[k + j * K];
            C[i + j * M] = total;
        }
    }
#elif defined(__AVX2__)
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            __m256d sum = _mm256_setzero_pd();
            size_t k = 0;
            for (; k + 3 < K; k += 4) {
                __m256d a = _mm256_loadu_pd(A + i + k * M);
                __m256d b = _mm256_loadu_pd(B + k + j * K);
                sum = _mm256_fmadd_pd(a, b, sum);
            }
            double buffer[4];
            _mm256_storeu_pd(buffer, sum);
            double total = 0.0;
            for (int z = 0; z < 4; ++z) total += buffer[z];
            for (; k < K; ++k)
                total += A[i + k * M] * B[k + j * K];
            C[i + j * M] = total;
        }
    }
#else
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N; ++j)
        for (size_t i = 0; i < M; ++i) {
            double sum = 0;
            for (size_t k = 0; k < K; ++k)
                sum += A[i + k * M] * B[k + j * K];
            C[i + j * M] = sum;
        }
#endif
}

template <typename T>
const Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) {
    if (this->sz_cols != B.sz_rows)
        throw std::runtime_error("Dimension mismatch for matrix multiplication");

    size_t M = this->sz_rows;
    size_t N = B.sz_cols;
    size_t K = this->sz_cols;

    Matrix<T> C(M, N);

    if constexpr (std::is_same<T, float>::value) {
        const char trans = 'N';
        int m = static_cast<int>(M), n = static_cast<int>(N), k = static_cast<int>(K);
        int lda = m, ldb = k, ldc = m;
        float alpha = 1.0f, beta = 0.0f;
        sgemm_(&trans, &trans, &m, &n, &k,
               &alpha, this->data, &lda, B.data, &ldb, &beta, C.data, &ldc);
    } else if constexpr (std::is_same<T, double>::value) {
        const char trans = 'N';
        int m = static_cast<int>(M), n = static_cast<int>(N), k = static_cast<int>(K);
        int lda = m, ldb = k, ldc = m;
        double alpha = 1.0, beta = 0.0;
        dgemm_(&trans, &trans, &m, &n, &k,
               &alpha, this->data, &lda, B.data, &ldb, &beta, C.data, &ldc);
    } else {
        matmul_simd_fallback(this->data, B.data, C.data, M, N, K);
    }

    return C;
}
    */