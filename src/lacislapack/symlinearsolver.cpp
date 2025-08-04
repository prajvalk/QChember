/*
 * NewSCF
 * Copyright (C) 2025, Prajval K
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "newscf/cis.hpp"
#include "lacislapack.hpp"

#include <cassert>

extern "C" {
        void dsysv_ (char* uplo, int* n, int* nrhs, double* a, int* lda,
                                           int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

        void dsysv_rook_ (char* uplo, int* n, int* nrhs, double* a, int* lda,
                                           int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

        void dsysv_rk_ (char* uplo, int* n, int* nrhs, double* a, int* lda,
                double* e, int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

        void dsysv_aa_ (char* uplo, int* n, int* nrhs, double* a, int* lda,
                                           int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

        void dsysv_aa_2stage_ (char* uplo, int* n, int* nrhs, double* a, int* lda, double* tb, int* ltb,
                                           int* ipiv, int* ipiv2, double* b, int* ldb, double* work, int* lwork, int* info);
}

namespace newscf::cis {

#define _IMPLEMENT_LINEAR_SOLVER_SET_UPLO(TYPE) \
        template<> \
        void LinearSolver<TYPE>::set_uplo(const char ch) { \
                CH1 = ch; \
        }

#define _IMPLEMENT_LINEAR_SOLVER_SET_N(TYPE) \
        template<> \
        void LinearSolver<TYPE>::set_n (const int n) { \
                I1 = n; \
                I3 = n; \
                I4 = n; \
        }

#define _IMPLEMENT_LINEAR_SOLVER_SET_NRHS(TYPE) \
        template<> \
        void LinearSolver<TYPE>::set_nrhs (const int nrhs) { \
                I2 = nrhs; \
        }

#define _IMPLEMENT_LINEAR_SOLVER_SET_LDA(TYPE) \
        template<> \
        void LinearSolver<TYPE>::set_lda (const int lda) { \
                I3 = lda; \
        }

#define _IMPLEMENT_LINEAR_SOLVER_SET_LDB(TYPE) \
        template<> \
        void LinearSolver<TYPE>::set_ldb (const int ldb) { \
                I4 = ldb; \
        }

#define _LAPACK_BACKEND_IMPLEMENT_LINEAR_SOLVER(TYPE) \
        _IMPLEMENT_LINEAR_SOLVER_SET_UPLO(TYPE) \
        _IMPLEMENT_LINEAR_SOLVER_SET_N(TYPE) \
        _IMPLEMENT_LINEAR_SOLVER_SET_NRHS(TYPE) \
        _IMPLEMENT_LINEAR_SOLVER_SET_LDA(TYPE) \
        _IMPLEMENT_LINEAR_SOLVER_SET_LDB(TYPE)

#ifndef LACIS_DISABLE_CHECKS
#define _CHECK_LINEAR_SOLVER_PARAMETERS(I1, I2, I3, I4) \
        assert(I1 > 0); \
        assert(I2 >= 1); \
        assert(I3 >= I1); \
        assert(I4 >= I1);
#elif
#define _CHECK_LINEAR_SOLVER_PARAMETERS(I1, I2, I3, I4)
#endif

#define _IMPLEMENT_DSYSV_LIKE_INIT_LINEAR_SOLVER(TYPE, LAPACKSYMB) \
        template<> \
        void LinearSolver<TYPE>::init() { \
                size_t workbytes = 0; \
                work_ptrs = MALLOC(sizeof(size_t) * 4); \
                _CHECK_LINEAR_SOLVER_PARAMETERS(I1, I2, I3, I4) \
                LWORK = -1; \
                work = MALLOC(sizeof(double) * 1); \
                LAPACKSYMB (&CH1, &I1, &I2, nullptr, &I3, nullptr, nullptr, &I4, reinterpret_cast<double*>(work), &LWORK, &INFO); \
                _LAPACK_INFO_CHECK(INFO) \
                LWORK = static_cast<int>(reinterpret_cast<double*>(work)[0]); \
                FREE(work); \
                reinterpret_cast<size_t*>(work_ptrs)[0] = workbytes; \
                workbytes += sizeof(double) * I3 * I1; \
                reinterpret_cast<size_t*>(work_ptrs)[1] = workbytes; \
                workbytes += sizeof(double) * I4 * I2; \
                reinterpret_cast<size_t*>(work_ptrs)[2] = workbytes; \
                workbytes += sizeof(int)    * I1; \
                reinterpret_cast<size_t*>(work_ptrs)[3] = workbytes; \
                workbytes += sizeof(double) * LWORK; \
                work = MALLOC(workbytes); \
        }

        template<>
        void LinearSolver<DSYSV_RK>::init() {
                size_t workbytes = 0;
                work_ptrs = MALLOC(sizeof(size_t) * 5);
                _CHECK_LINEAR_SOLVER_PARAMETERS(I1, I2, I3, I4)
                LWORK = -1;
                work = MALLOC(sizeof(double) * 1);
                dsysv_rk_ (&CH1, &I1, &I2, nullptr, &I3, nullptr, nullptr, nullptr, &I4, reinterpret_cast<double*>(work), &LWORK, &INFO); \
                _LAPACK_INFO_CHECK(INFO) \
                LWORK = static_cast<int>(reinterpret_cast<double*>(work)[0]);
                FREE(work);
                reinterpret_cast<size_t*>(work_ptrs)[0] = workbytes;
                workbytes += sizeof(double) * I3 * I1;
                reinterpret_cast<size_t*>(work_ptrs)[1] = workbytes;
                workbytes += sizeof(double) * I4 * I2;
                reinterpret_cast<size_t*>(work_ptrs)[2] = workbytes;
                workbytes += sizeof(int)    * I1;
                reinterpret_cast<size_t*>(work_ptrs)[3] = workbytes;
                workbytes += sizeof(double) * LWORK;
                reinterpret_cast<size_t*>(work_ptrs)[4] = workbytes;
                workbytes += sizeof(double) * I1;
                work = MALLOC(workbytes);
        }

        template<>
        void LinearSolver<DSYSV_AA_2>::init() {
                size_t workbytes = 0;
                work_ptrs = MALLOC(sizeof(size_t) * 6);
                _CHECK_LINEAR_SOLVER_PARAMETERS(I1, I2, I3, I4)
                LWORK = -1;
                I5 = -1;
                work = MALLOC(sizeof(double) * 2);
                double* tb = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + sizeof(double));
                dsysv_aa_2stage_ (&CH1, &I1, &I2, nullptr, &I3, tb, &I5, nullptr, nullptr, nullptr, &I4, reinterpret_cast<double*>(work), &LWORK, &INFO); \
                _LAPACK_INFO_CHECK(INFO)
                LWORK = static_cast<int>(reinterpret_cast<double*>(work)[0]);
                I5 = static_cast<int>(reinterpret_cast<double*>(work)[1]);
                FREE(work);
                reinterpret_cast<size_t*>(work_ptrs)[0] = workbytes;
                workbytes += sizeof(double) * I3 * I1;
                reinterpret_cast<size_t*>(work_ptrs)[1] = workbytes;
                workbytes += sizeof(double) * I4 * I2;
                reinterpret_cast<size_t*>(work_ptrs)[2] = workbytes;
                workbytes += sizeof(int)    * I1;
                reinterpret_cast<size_t*>(work_ptrs)[3] = workbytes;
                workbytes += sizeof(double) * LWORK;
                reinterpret_cast<size_t*>(work_ptrs)[4] = workbytes; // TB
                workbytes += sizeof(double) * I5;
                reinterpret_cast<size_t*>(work_ptrs)[5] = workbytes; // IPIV2
                workbytes += sizeof(int) * I1;
                work = MALLOC(workbytes);
        }

#define _EXPOSE_LACIS_LINEAR_SOLVER_LAPACK_INTERFACE(T, E, R) \
        template<> \
        template<> \
        void LinearSolver<T>::set_a<E, R> (const LACIS<E, R>& A) { \
                const void* cA = A.data; \
                void* workA = (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]); \
                MEMCPY(workA, cA, sizeof(R) * I3 * I1); \
        } \
        template<> \
        template<> \
        void LinearSolver<T>::set_b<E, R> (const LACIS<E, R>& B) { \
                const void* cB = B.data; \
                void* workB = (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]); \
                MEMCPY(workB, cB, sizeof(R) * I4 * I2); \
        } \
        template<> \
        template<> \
        void LinearSolver<T>::copy_result (LACIS<E, R>& result) { \
                void* res = reinterpret_cast<void*>(result.data); \
                const void* workB = (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]); \
                MEMCPY(res, workB, sizeof(R) * I2 * I4); \
        }

#define _IMPLEMENT_LACIS_LAPACK_SOLVER(TYPE, LAPACKSYMB) \
        template<> \
        void LinearSolver<TYPE>::solve() { \
                double* workA     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]); \
                double* workB     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]); \
                int*    workIPIV  = reinterpret_cast<int*>   (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[2]); \
                double* workspace = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[3]); \
                LAPACKSYMB (&CH1, &I1, &I2, workA, &I3, workIPIV, workB, &I4, workspace, &LWORK, &INFO); \
                _LAPACK_INFO_CHECK(INFO) \
        }

        template<>
        void LinearSolver<DSYSV_RK>::solve() {
                double* workA     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]);
                double* workB     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]);
                int*    workIPIV  = reinterpret_cast<int*>   (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[2]);
                double* workspace = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[3]);
                double* workE     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[4]);
                dsysv_rk_ (&CH1, &I1, &I2, workA, &I3, workE, workIPIV, workB, &I4, workspace, &LWORK, &INFO); \
                _LAPACK_INFO_CHECK(INFO) \
        }

        template<>
        void LinearSolver<DSYSV_AA_2>::solve() {
                double* workA     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]);
                double* workB     = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]);
                int*    workIPIV  = reinterpret_cast<int*>   (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[2]);
                double* workspace = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[3]);
                double* workTB    = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[4]);
                int*    workIPIV2 = reinterpret_cast<int*>   (reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[5]);
                dsysv_aa_2stage_ (&CH1, &I1, &I2, workA, &I3, workTB, &I5, workIPIV, workIPIV2, workB, &I4, workspace, &LWORK, &INFO);
                _LAPACK_INFO_CHECK(INFO) \
        }

#define _IMPLEMENT_ALL_LACIS_LAPACK_DSYSV_LIKE_INERFACE(SOLTYPE, LACIS_TYPE, LAPACKSYMB, DT) \
        _LAPACK_BACKEND_IMPLEMENT_LINEAR_SOLVER(SOLTYPE) \
        _IMPLEMENT_DSYSV_LIKE_INIT_LINEAR_SOLVER(SOLTYPE, LAPACKSYMB) \
        _EXPOSE_LACIS_LINEAR_SOLVER_LAPACK_INTERFACE(SOLTYPE, LACIS_TYPE, DT) \
        _IMPLEMENT_LACIS_LAPACK_SOLVER(SOLTYPE, LAPACKSYMB)

#define _IMPLEMENT_LIMITED_LACIS_LAPACK_DSYSV_LIKE_INTERFACE(SOLTYPE, LACIS_TYPE, DT) \
        _LAPACK_BACKEND_IMPLEMENT_LINEAR_SOLVER(SOLTYPE) \
        _EXPOSE_LACIS_LINEAR_SOLVER_LAPACK_INTERFACE(SOLTYPE, LACIS_TYPE, DT)

        _IMPLEMENT_ALL_LACIS_LAPACK_DSYSV_LIKE_INERFACE(DSYSV, DENSE_LAPACK, dsysv_, double)
        _IMPLEMENT_ALL_LACIS_LAPACK_DSYSV_LIKE_INERFACE(DSYSV_ROOK, DENSE_LAPACK, dsysv_rook_, double)
        _IMPLEMENT_LIMITED_LACIS_LAPACK_DSYSV_LIKE_INTERFACE(DSYSV_RK, DENSE_LAPACK, double)
        _IMPLEMENT_ALL_LACIS_LAPACK_DSYSV_LIKE_INERFACE(DSYSV_AA, DENSE_LAPACK, dsysv_aa_, double)
        _IMPLEMENT_LIMITED_LACIS_LAPACK_DSYSV_LIKE_INTERFACE(DSYSV_AA_2, DENSE_LAPACK, double)

}