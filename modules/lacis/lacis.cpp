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

#include "qchember/api/cis.hpp"

#include <cassert>

namespace newscf::cis {

    template<LACIS_TYPE BACKEND, typename T>
    LACIS<BACKEND, T>::LACIS() {
        data = nullptr;
        dims = nullptr;
        ndims = 0;
        ndata = 0;
    }

    template<LACIS_TYPE BACKEND, typename T>
    LACIS<BACKEND, T>::~LACIS() {
        FREE(data);
        FREE(dims);
    }

    template<LACIS_TYPE BACKEND, typename T>
    LACIS<BACKEND, T>& LACIS<BACKEND, T>::operator=(const LACIS<BACKEND, T>& rhs) {
        if (this != &rhs) {
            FREE(data);
            FREE(dims);
            ndata = rhs.ndata;
            ndims = rhs.ndims;
            data = MALLOC(sizeof(double) * ndata);
            MEMCPY(data, rhs.data, sizeof(double) * ndata);
            dims = MALLOC(sizeof(size_t) * ndims);
            MEMCPY(dims, rhs.dims, sizeof(size_t) * ndims);
        }
        return *this;
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToVector(const size_t size) {
        dims = MALLOC(sizeof(size_t) * 1);
        data = MALLOC(sizeof(T) * size);
        reinterpret_cast<size_t*>(dims)[0] = size;
        ndata = size;
        ndims = 1;
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToMatrix(const size_t r, const size_t c) {
        dims = MALLOC(sizeof(size_t) * 2);
        data = MALLOC(sizeof(T) * r * c);
        reinterpret_cast<size_t*>(dims)[0] = r;
        reinterpret_cast<size_t*>(dims)[1] = c;
        ndata = r * c;
        ndims = 2;
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToMatrix(const size_t i) {
        resizeToMatrix(i, i);
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToTensor3D(const size_t i, const size_t j, const size_t k) {
        dims = MALLOC(sizeof(size_t) * 3);
        data = MALLOC(sizeof(T) * i * j * k);
        reinterpret_cast<size_t*>(dims)[0] = i;
        reinterpret_cast<size_t*>(dims)[1] = j;
        reinterpret_cast<size_t*>(dims)[2] = k;
        ndata = i * j * k;
        ndims = 3;
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToTensor3D(const size_t i) {
        resizeToTensor3D(i, i, i);
    }


    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToTensor4D(const size_t i, const size_t j, const size_t k, const size_t l) {
        dims = MALLOC(sizeof(size_t) * 4);
        data = MALLOC(sizeof(T) * i * j * k * l);
        reinterpret_cast<size_t*>(dims)[0] = i;
        reinterpret_cast<size_t*>(dims)[1] = j;
        reinterpret_cast<size_t*>(dims)[2] = k;
        reinterpret_cast<size_t*>(dims)[3] = l;
        ndata = i * j * k * l;
        ndims = 4;
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::resizeToTensor4D(const size_t i) {
        resizeToTensor4D(i, i, i, i);
    }


    template<LACIS_TYPE BACKEND, typename T>
    T LACIS<BACKEND, T>::vectorGet (const size_t i) {
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(i < ndata);
        assert(ndims >= 1);
#endif
        return reinterpret_cast<T*>(data)[i];
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::vectorSet (const size_t i, const T val) {
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(i < ndata);
        assert(ndims >= 1);
#endif
        reinterpret_cast<T*>(data)[i] = val;
    }

    template<LACIS_TYPE BACKEND, typename T>
    T LACIS<BACKEND, T>::matrixGet (const size_t i, const size_t j) {
        const size_t rows = reinterpret_cast<size_t*>(dims)[0];
        const size_t cols = reinterpret_cast<size_t*>(dims)[1];
        const size_t idx  = j * rows + i;
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(idx < ndata);
        assert(ndims >= 2);
        assert(i < rows);
        assert(j < cols);
#endif
        return reinterpret_cast<T*>(data)[idx];
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::matrixSet (const size_t i, const size_t j, const T val) {
        const size_t rows = reinterpret_cast<size_t*>(dims)[0];
        const size_t cols = reinterpret_cast<size_t*>(dims)[1];
        const size_t idx  = j * rows + i;
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(idx < ndata);
        assert(ndims >= 2);
        assert(i < rows);
        assert(j < cols);
#endif
        reinterpret_cast<T*>(data)[idx] = val;
    }

    template<LACIS_TYPE BACKEND, typename T>
    T LACIS<BACKEND, T>::tensor3DGet (const size_t i, const size_t j, const size_t k) {
        const size_t ld1 = reinterpret_cast<size_t*>(dims)[0];
        const size_t ld2 = reinterpret_cast<size_t*>(dims)[1];
        const size_t ld3 = reinterpret_cast<size_t*>(dims)[2];
        const size_t idx = k * ld2 * ld1 + j * ld1 + i;
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(idx < ndata);
        assert(ndims >= 3);
        assert(i < ld1);
        assert(j < ld2);
        assert(k < ld3);
#endif
        return reinterpret_cast<T*>(data)[idx];
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::tensor3DSet (const size_t i, const size_t j, const size_t k, const T val) {
        const size_t ld1 = reinterpret_cast<size_t*>(dims)[0];
        const size_t ld2 = reinterpret_cast<size_t*>(dims)[1];
        const size_t ld3 = reinterpret_cast<size_t*>(dims)[2];
        const size_t idx = k * ld2 * ld1 + j * ld1 + i;
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(idx < ndata);
        assert(ndims >= 3);
        assert(i < ld1);
        assert(j < ld2);
        assert(k < ld3);
#endif
        reinterpret_cast<T*>(data)[idx] = val;
    }

    template<LACIS_TYPE BACKEND, typename T>
    T LACIS<BACKEND, T>::tensor4DGet (const size_t i, const size_t j, const size_t k, const size_t l) {
        const size_t ld1 = reinterpret_cast<size_t*>(dims)[0];
        const size_t ld2 = reinterpret_cast<size_t*>(dims)[1];
        const size_t ld3 = reinterpret_cast<size_t*>(dims)[2];
        const size_t ld4 = reinterpret_cast<size_t*>(dims)[3];
        const size_t idx = l * ld3 * ld2 * ld1 + k * ld2 * ld1 + j * ld1 + i;
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(idx < ndata);
        assert(ndims >= 4);
        assert(i < ld1);
        assert(j < ld2);
        assert(k < ld3);
        assert(l < ld4);
#endif
        return reinterpret_cast<T*>(data)[idx];
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::tensor4DSet (const size_t i, const size_t j, const size_t k, const size_t l, const T val) {
        const size_t ld1 = reinterpret_cast<size_t*>(dims)[0];
        const size_t ld2 = reinterpret_cast<size_t*>(dims)[1];
        const size_t ld3 = reinterpret_cast<size_t*>(dims)[2];
        const size_t ld4 = reinterpret_cast<size_t*>(dims)[3];
        const size_t idx = l * ld3 * ld2 * ld1 + k * ld2 * ld1 + j * ld1 + i;
#ifndef LACIS_DISABLE_BOUNDS_CHECKS
        assert(idx < ndata);
        assert(ndims >= 4);
        assert(i < ld1);
        assert(j < ld2);
        assert(k < ld3);
        assert(l < ld4);
#endif
        reinterpret_cast<T*>(data)[idx] = val;
    }

    template <LACIS_TYPE BACKEND, typename T>
    template <NormMetric type>
    T LACIS<BACKEND, T>::compareTo(const LACIS<BACKEND, T>& lhs) {
        T res = T(0);
        if (this->ndata != lhs.ndata)
            HANDLE_ERROR("Incompatible comparison", 90);

        if constexpr (type == L2) {
            for (size_t i = 0; i < this->ndata; i++)
                res += pow(this->vectorGet(i) - lhs.vectorGet(i), 2);
            return sqrt(res);
        } else if constexpr (type == L1) {
            for (size_t i = 0; i < this->ndata; i++)
                res += std::abs(this->vectorGet(i) - lhs.vectorGet(i));
            return res;
        } else {
            HANDLE_ERROR("Unsupported norm type", 91);
        }
        return res;
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::fill(T t) {
        MEMSET(data, static_cast<int>(t), ndata * sizeof(T));
    }

    template<LACIS_TYPE BACKEND, typename T>
    void LACIS<BACKEND, T>::vectorPrint() {
        for (size_t i = 0; i < this->ndata; i++)
            if (this->vectorGet(i) != 0)
                std::cout << i << " " << this->vectorGet(i) << "\n";
    }

    template class LACIS<DENSE_LAPACK, double>;

}