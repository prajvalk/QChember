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
    void dsygv_(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
                 int* info);
}

namespace newscf::cis {

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_itype(const int itype) {
        I1 = itype;
    }

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_jobz(const char jobz) {
        CH1 = jobz;
    }

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_uplo(const char uplo) {
        CH2 = uplo;
    }

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_n(const int n) {
        I2 = n;
        I3 = n;
        I4 = n;
    }

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_lda(const int lda) {
        I2 = lda;
    }

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_ldb(const int ldb) {
        I3 = ldb;
    }

#define _CHECK_EIGENSOLVER_PARAMETERS(I1, I2, I3) \
    assert(I1 >= 1); \
    assert(I2 >= I1); \
    assert(I3 >= I1);

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::init() {
        size_t workbytes = 0;
        _CHECK_EIGENSOLVER_PARAMETERS(I1, I2, I3)
        work_ptrs = MALLOC(sizeof(size_t) * 5);
        LWORK1 = -1;
        INFO = 0;
        work = MALLOC(sizeof(double) * 1);
        dsygv_ (&I1, &CH1, &CH2, &I2, nullptr, &I3, nullptr, &I4, nullptr, reinterpret_cast<double*>(work), &LWORK1, &INFO);
        _LAPACK_INFO_CHECK(INFO);
        LWORK1 = static_cast<int>(reinterpret_cast<double*>(work)[0]);
        FREE(work);
        reinterpret_cast<size_t*>(work_ptrs)[0] = workbytes;
        workbytes += sizeof(double) * I3 * I2;
        reinterpret_cast<size_t*>(work_ptrs)[1] = workbytes;
        workbytes += sizeof(double) * I4 * I2;
        reinterpret_cast<size_t*>(work_ptrs)[2] = workbytes;
        workbytes += sizeof(double) * I2;
        reinterpret_cast<size_t*>(work_ptrs)[3] = workbytes;
        workbytes += sizeof(double) * LWORK1;
        work = MALLOC(workbytes);
    }

    template<>
    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_a(const LACIS<DENSE_LAPACK, double> &A) {
        const void* cA = A.data;
        void* workA = reinterpret_cast<void*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]);
        MEMCPY(workA, cA, sizeof(double) * I2 * I3);
    }

    template<>
    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::set_b(const LACIS<DENSE_LAPACK, double> &B) {
        const void* cB = B.data;
        void* workB = reinterpret_cast<void*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]);
        MEMCPY(workB, cB, sizeof(double) * I2 * I4);
    }

    template<>
    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::copy_eigenvalues(LACIS<DENSE_LAPACK, double> &res) {
        void* cres = res.data;
        const void* workEVAL = reinterpret_cast<void*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[2]);
        MEMCPY(cres, workEVAL, sizeof(double) * I2);
    }

    template<>
    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::copy_eigenvectors(LACIS<DENSE_LAPACK, double> &res) {
        void* cA = res.data;
        const void* workA = reinterpret_cast<void*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]);
        MEMCPY(cA, workA, sizeof(double) * I2 * I3);
    }

    template<>
    void GeneralizedEigenvalueSolver<DSYGV>::solve() {
        double* wA = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[0]);
        double* wB = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[1]);
        double* wE = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[2]);
        double* wW = reinterpret_cast<double*>(reinterpret_cast<char*>(work) + reinterpret_cast<size_t*>(work_ptrs)[3]);
        dsygv_ (&I1, &CH1, &CH2, &I2, wA, &I3, wB, &I4, wE, wW, &LWORK1, &INFO);
    }

}