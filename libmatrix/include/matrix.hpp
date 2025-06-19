#ifndef NEWSCF_LIBMATRIX_MATRIX_HPP
#define NEWSCF_LIBMATRIX_MATRIX_HPP

#include "commons.hpp"
    
#include <iostream>
#include <algorithm>
#include <cmath>

namespace newscf::libmatrix {

    // OpenMP+SIMD Intrinsics

    void add_float_simd            (const float* A, const float* B, float* C, size_t N);
    void add_inplace_float_simd    (float* C, const float* A, size_t N);
    void add_double_simd           (const double* A, const double* B, double* C, size_t N);
    void add_inplace_double_simd   (double* C, const double* A, size_t N);

    void sub_float_simd            (const float* A, const float* B, float* C, size_t N);
    void sub_inplace_float_simd    (float* C, const float* A, size_t N);
    void sub_double_simd           (const double* A, const double* B, double* C, size_t N);
    void sub_inplace_double_simd   (double* C, const double* A, size_t N);

    void mul_scalar_float_simd     (const float* A, float* C, float scalar, size_t N);
    void mul_scalar_inplace_float_simd (float* C, float scalar, size_t N);
    void mul_scalar_double_simd    (const double* A, double* C, double scalar, size_t N);
    void mul_scalar_inplace_double_simd (double* C, double scalar, size_t N);
    
    // OpenMP Generics

    template<typename T>
    inline void add_generic (const T* A, const T* B, T* C, size_t N) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] = A[i] + B[i];
    }

    template<typename T>
    inline void add_inplace_generic (T* A, const T* B, size_t N) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            A[i] += B[i];
    }

    template<typename T>
    inline void sub_generic (const T* A, const T* B, T* C, size_t N) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] = A[i] - B[i];
    }

    template<typename T>
    inline void sub_inplace_generic (T* A, const T* B, size_t N) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            A[i] -= B[i];
    }

    template <typename T>
    inline void mul_scalar_generic(const T* A, T* C, const T scalar, size_t N) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] = A[i] * scalar;
    }

    template <typename T>
    inline void mul_scalar_inplace_generic(T* C, const T scalar, size_t N) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i)
            C[i] *= scalar;
    }

    template <typename T>
    class Matrix {
    public:
        unsigned int sz_rows; // number of rows
        unsigned int sz_cols; // number of cols
        T* data;              // martrix structure, stored in column-major format
        Matrix() {
            data = new T[1];
            sz_rows = 1;
            sz_cols = 1;
        }

        ~Matrix() {
            delete[] data;
            sz_rows = 0;
            sz_cols = 0;
        }

        Matrix(const Matrix& oth) {
            data = new T [oth.sz_rows * oth.sz_cols];
            std::copy(oth.data, oth.data + oth.sz_cols * oth.sz_rows, data);
            sz_rows = oth.sz_rows;
            sz_cols = oth.sz_cols;
        }

        Matrix (const unsigned int r, const unsigned int c) {
            data = new T [r * c];
            sz_rows = r;
            sz_cols = c;
        }

        explicit Matrix(T* ref, const unsigned int r, const unsigned int c) {
            data = new T[r * c];
            std::copy(ref, ref + r * c, data);
            sz_rows = r;
            sz_cols = c;
        }

        Matrix& operator = (const Matrix& oth) {
            if (this != &oth) {
                T* cpy = new T [oth.sz_cols * oth.sz_rows];
                std::copy (oth.data, oth.data + oth.sz_cols * oth.sz_rows, cpy);
                delete[] data;
                data = cpy;
                sz_rows = oth.sz_rows;
                sz_cols = oth.sz_cols;
            }
            return *this;
        }

        const Matrix<T>  operator+  (const Matrix& other) {
            if (this->sz_rows != other.sz_rows || this->sz_cols != other.sz_cols)
                throw std::runtime_error("Matrix dimensions do not match");

            size_t M = this->sz_rows;
            size_t N = this->sz_cols;
            size_t total = M * N;

            Matrix<T> result(M, N);

            const T* A_ptr = this->data;
            const T* B_ptr = other.data;
            T* C_ptr = result.data;

            if constexpr (std::is_same<T, float>::value) {
                add_float_simd (A_ptr, B_ptr, C_ptr, total);
            } else if constexpr (std::is_same<T, double>::value) {
                add_double_simd (A_ptr, B_ptr, C_ptr, total);
            } else {
                add_generic (A_ptr, B_ptr, C_ptr, total);
            }

            return result;
        }

        const Matrix<T> &operator+= (const Matrix<T> &other) {
            if (this->sz_rows != other.sz_rows || this->sz_cols != other.sz_cols)
                throw std::runtime_error("Matrix dimensions do not match");
            size_t M = this->sz_rows;
            size_t N = this->sz_cols;
            size_t total = M * N;

            T* C_ptr = this->data;
            T* A_ptr = other.data;
            if constexpr (std::is_same<T, float>::value) {
                add_inplace_float_simd (C_ptr, A_ptr, total);
            } else if constexpr (std::is_same<T, double>::value) {
                add_inplace_double_simd (C_ptr, A_ptr, total);
            } else {
                add_inplace_generic (C_ptr, A_ptr, total);
            }

            return *this;
        }

        const Matrix<T>  operator-  (const Matrix& other) {
            if (this->sz_rows != other.sz_rows || this->sz_cols != other.sz_cols)
                throw std::runtime_error("Matrix dimensions do not match");

            size_t M = this->sz_rows;
            size_t N = this->sz_cols;
            size_t total = M * N;

            Matrix<T> result(M, N);

            const T* A_ptr = this->data;
            const T* B_ptr = other.data;
            T* C_ptr = result.data;

            if constexpr (std::is_same<T, float>::value) {
                sub_float_simd (A_ptr, B_ptr, C_ptr, total);
            } else if constexpr (std::is_same<T, double>::value) {
                sub_double_simd (A_ptr, B_ptr, C_ptr, total);
            } else {
                sub_generic (A_ptr, B_ptr, C_ptr, total);
            }

            return result;
        }

        const Matrix<T>  operator-= (const Matrix<T> &other) {
            if (this->sz_rows != other.sz_rows || this->sz_cols != other.sz_cols)
                throw std::runtime_error("Matrix dimensions do not match");

            size_t M = this->sz_rows;
            size_t N = this->sz_cols;
            size_t total = M * N;
            T* C_ptr = this->data;
            T* A_ptr = other.data;
            if constexpr (std::is_same<T, float>::value) {
                sub_inplace_float_simd (C_ptr, A_ptr, total);
            } else if constexpr (std::is_same<T, double>::value) {
                sub_inplace_double_simd (C_ptr, A_ptr, total);
            } else {
                sub_inplace_generic (C_ptr, A_ptr, total);
            }

            return *this;
        }

        const Matrix<T> operator*  (const T s) {
            size_t M = this->sz_rows;
            size_t N = this->sz_cols;
            size_t total = M * N;

            Matrix<T> result(M, N);

            const T* A_ptr = this->data;
            T* C_ptr = result.data;

            if constexpr (std::is_same<T, float>::value) {
                mul_scalar_float_simd (A_ptr, C_ptr, s, total);
            } else if constexpr (std::is_same<T, double>::value) {
                mul_scalar_double_simd (A_ptr, C_ptr, s, total);
            } else {
                mul_scalar_generic (A_ptr, C_ptr, s, total);
            }

            return result;
        }
        const Matrix<T> operator*= (const T s) {
            size_t M = this->sz_rows;
            size_t N = this->sz_cols;
            size_t total = M * N;

            T* A_ptr = this->data;

            if constexpr (std::is_same<T, float>::value) {
                mul_scalar_inplace_float_simd (A_ptr, s, total);
            } else if constexpr (std::is_same<T, double>::value) {
                mul_scalar_inplace_double_simd (A_ptr, s, total);
            } else {
                mul_scalar_inplace_generic (A_ptr, s, total);
            }

            return *this;
        }

        /* {
                    if (sz_rows != A.sz_rows || sz_cols != A.sz_cols) {
                        std::cerr << "  [ERR] Invalid shape" << std::endl;
                        exit(-1);
                    }
                    T* ndat = new T[sz_cols * sz_rows];
                    for(auto i = 0; i < sz_cols * sz_rows; i++) ndat[i] = data[i] + A.data[i];
                    Matrix res (ndat, sz_rows, sz_cols);
                    delete[] ndat;
                    return res;
                }*/

        /*Matrix operator - (const Matrix& A) const {
            if (sz_rows != A.sz_rows || sz_cols != A.sz_cols) {
                std::cerr << "  [ERR] Invalid shape" << std::endl;
                exit(-1);
            }
            T* ndat = new T[sz_cols * sz_rows];
            for(auto i = 0; i < sz_cols * sz_rows; i++) ndat[i] = data[i] - A.data[i];
            Matrix res (ndat, sz_rows, sz_cols);
            delete[] ndat;
            return res;
        }*/

        /*template <typename E>
        Matrix operator * (E e) const {
            T* ndat = new T[sz_cols * sz_rows];
            for(auto i = 0; i < sz_cols * sz_rows; i++) ndat[i] = data[i] * e;
            Matrix res (ndat, sz_rows, sz_cols);
            delete[] ndat;
            return res;
        }*/

        inline void set(const unsigned int i, const unsigned int j, const T t) const {
            if (i >= sz_rows || j >= sz_cols) {
                std::cerr << "  [ERR] : Invalid Memory Access \n";
                std::cerr << "  [ERR] : Valid Bounds ("<< sz_rows << "," << sz_cols << ") \n";
                std::cerr << "  [ERR] : Accessed ("<< i << "," << j << ") \n";
            } else
                data[j * sz_rows + i] = t;
        }

        inline void clean() const {
            for (unsigned int i = 0; i < sz_rows; i++)
                for (unsigned int j = 0; j < sz_cols; j++)
                    data[j * sz_rows + i] = 0;
        }

        inline T get(const unsigned int i, const unsigned int j) const {
            if (i >= sz_rows || j >= sz_cols) {
                std::cerr << "  [ERR] : Invalid Memory Access \n";
                std::cerr << "  [ERR] : Valid Bounds ("<< sz_rows << "," << sz_cols << ") \n";
                std::cerr << "  [ERR] : Accessed ("<< i << "," << j << ") \n";
                return -1;
            } else
                return data[j * sz_rows + i];
        }

        inline void transpose() {
            Matrix trans (sz_cols, sz_rows);
            for (auto i = 0; i < sz_rows; i++)
                for (auto j = 0; j < sz_cols; j++)
                    trans.set(j, i, get(i, j));
            (*this) = trans;
        }

        inline unsigned int countNonZero() {
            unsigned int counter = 0;
            for (unsigned int i = 0; i < sz_rows * sz_cols; i++) if (data[i] != 0) counter++;
            return counter;
        }

        inline void print() const {
            for (unsigned int i = 0; i < sz_rows; i++) {
                for (unsigned int j = 0; j < sz_cols; j++)
                    std::cout << data[j * sz_rows + i] << "\t";
                std::cout << "\n \n";
            }
        }

        inline double compareTo(const Matrix& m, const NORM err_method) const {
            if (sz_rows == m.sz_rows && sz_cols == m.sz_cols) {
                if (err_method == LInf) {
                    double max_err = 0;
                    for (auto i = 0; i < sz_cols * sz_rows; i++) {
                        double err = fabs(data[i] - m.data[i]);
                        if (err > max_err) max_err = err;
                    }
                    return max_err;
                } else if (err_method == L2) {
                    double sse = 0;
                    for (auto i = 0; i < sz_cols * sz_rows; i++) {
                        double err = fabs(data[i] - m.data[i]);
                        sse += err * err;
                    }
                    return sqrt(sse);
                } else if (err_method == L1) {
                    double sae = 0;
                    for (auto i = 0; i < sz_cols * sz_rows; i++) {
                        double err = fabs(data[i] - m.data[i]);
                        sae += err;
                    }
                    return sae;
                } else if (err_method == L0) {
                    double sum = 0;
                    for (auto i = 0; i < sz_cols * sz_rows; i++) {
                        double err = fabs(data[i] - m.data[i]);
                        sum += err;
                    }
                    return sum / (sz_cols * sz_rows);
                } /* 
                
                :(

                else if (err_method == MEDIAN) {
                    double* errs = new double[sz_cols * sz_rows];
                    double median = 0;
                    for (auto i = 0; i < sz_cols * sz_rows; i++) errs[i] = fabs(data[i] - m.data[i]);
                    std::sort(errs, errs + sz_cols * sz_rows);
                    if ((sz_cols * sz_rows) % 2 == 0)
                        median = errs[sz_cols * sz_rows / 2];
                    else
                        median = 0.5 * errs[sz_cols * sz_rows / 2] + 0.5 *errs[(sz_cols * sz_rows - 1) / 2];
                    delete[] errs;
                    return median;
                }*/
            } else {
                std::cerr << "  [ERR] : Comparision not possible between non-similar shape matrices. \n";
                return -1;
            }
            return -1;
        }
    };

} // namespace newscf::libmatrix




#endif // NEWSCF_LIBMATRIX_MATRIX_HPP