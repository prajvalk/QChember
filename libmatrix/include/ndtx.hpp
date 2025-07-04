#ifndef NDTX_HPP
#define NDTX_HPP
#include <bits/streambuf_iterator.h>

namespace newscf::ndtx {

    /* Experimental Lightweight NDTensor Routines */
    template <typename T>
    class NDTX {
    public:
        T* data;
        size_t  ndata;
        size_t* dims;
        size_t  ndims;

        NDTX () {
            data = nullptr;
            dims = nullptr;
            ndims = 0;
            ndata = 0;
        }

        ~NDTX () {
            delete[] data;
            delete[] dims;
        }

        NDTX& operator=(const NDTX<T>& other) {
            if (this != &other) {
                delete[] data;
                delete[] dims;
                ndims = other.ndims;
                ndata = other.ndata;
                dims  = new size_t[ndims];
                std::copy(other.dims, other.dims + ndims, dims);
                data  = new T[ndata];
                std::copy(other.data, other.data + ndata, data);
            }
            return *this;
        }

        inline void resizeToVector(const size_t size) {
            delete[] data;
            delete[] dims;
            dims = new size_t [1] {size};
            data = new T[size];
            ndata = size;
            ndims = 1;
        }

        inline void resizeToMatrix(const size_t rows, const size_t cols) {
            delete[] data;
            delete[] dims;
            dims = new size_t [2] {rows, cols};
            data = new T[rows * cols];
            ndata = rows * cols;
            ndims = 2;
        }

        inline void resizeToMatrix(const size_t size) {
            resizeToMatrix(size, size);
        }

        inline void resizeToTensor3D(const size_t i, const size_t j, const size_t k) {
            delete[] data;
            delete[] dims;
            dims = new size_t [3] {i, j, k};
            data = new T[i * j * k];
            ndata = i * j * k;
            ndims = 3;
        }

        inline void resizeToTensor3D(const size_t i) {
            resizeToTensor3D(i, i, i);
        }

        inline void resizeToTensor4D(const size_t i, const size_t j, const size_t k, const size_t l) {
            delete[] data;
            delete[] dims;
            dims = new size_t [4] {i, j, k, l};
            data = new T[i * j * k * l];
            ndata = i * j * k * l;
            ndims = 4;
        }

        inline void resizeToTensor4D(const size_t i) {
            resizeToTensor4D(i, i, i, i);
        }

        NDTX (const int ndims, size_t* dims) {
            this->dims = new size_t[ndims];
            std::copy(dims, dims + ndims, this->dims);
            this->ndims = ndims;
            this->ndata = 0;
            for (int i = 0; i < ndims; ++i) this->ndata *= dims[i];
            this->data = new T[this->ndata];
        }

        NDTX (const NDTX<T>& other) {
            data = new T [other.ndata];
            dims = new size_t[other.ndims];
            std::copy(other.dims, other.dims + other.ndims, dims);
            std::copy(other.data, other.data + other.ndata, data);
            this->ndata = other.ndata;
            this->ndims = other.ndims;
        }

        inline T vectorGet (const size_t idx) {
            return data[idx];
        }

        inline T matrixGet (const size_t i, const size_t j) {
            return data[j * dims[0] + i];
        }

        inline T tensor3DGet (const size_t i, const size_t j, const size_t k) {
            return data[k * dims[1] + j * dims[0] + i];
        }

        inline T tensor4DGet (const size_t i, const size_t j, const size_t k, const size_t l) {
            size_t index = ((i * dims[1] + j) * dims[2] + k) * dims[3] + l;
            return data[index];
        }

        inline void vectorSet (const size_t idx, T val) {
            data[idx] = val;
        }

        inline void matrixSet (const size_t i, const size_t j, const T val) {
            data[j * dims[0] + i] = val;
        }

        inline void tensor3DSet (const size_t i, const size_t j, const size_t k, const T val) {
            data[k * dims[1] + j * dims[0] + i] = val;
        }

        inline void tensor4DSet (const size_t i, const size_t j, const size_t k, const size_t l, const T val) {
            size_t index = ((i * dims[1] + j) * dims[2] + k) * dims[3] + l;
            data[index] = val;
        }

        inline size_t size () {
            return ndata;
        }

        inline size_t order () {
            return ndims;
        }

        inline void matrixDims (size_t& rows, size_t& cols) {
            rows = dims[0];
            cols = dims[1];
        }

        inline void tensor3DDims (size_t& i, size_t& j, size_t& k) {
            i = dims[0];
            j = dims[1];
            k = dims[2];
        }

        inline void tensor4DDims (size_t& i, size_t& j, size_t& k, size_t& l) {
            i = dims[0];
            j = dims[1];
            k = dims[2];
            l = dims[3];
        }

        inline void fill (const T val) {
            std::fill(data, data + ndata, val);
        }

        inline void fill (const T* vals, const size_t size) {
            std::copy (vals, vals + size, data);
        }

        inline void vectorPrint () {
            for (int i = 0; i < ndata; ++i)
                if (data[i] != 0)
                    std::cout << i << ", " << data[i] << std::endl;
        }

        inline void matrixPrint () {
            for (int i = 0; i < dims[0]; i++)
                for (int j = 0; j < dims[1]; j++)
                    if (this->matrixGet(i, j) != 0)
                        std::cout << i << ", " << j << ", " << this->matrixGet(i, j) << std::endl;
        }

    };

};

#endif
