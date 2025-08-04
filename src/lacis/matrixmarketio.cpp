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

#include <fstream>
#include <iomanip>
#include <vector>

namespace newscf::cis {

    template<LACIS_TYPE Type, typename T>
    int count_non_zero (LACIS<Type, T>& ref) {
        int count = 0;
        for (size_t i = 0; i < ref.ndata; ++i) {
            if (ref.vectorGet(i) != 0)
                ++count;
        }
        return count;
    }

    template<>
    template<>
    void LACIS<DENSE_LAPACK, double>::save_to_file<MATRIXMARKET> (const std::string& file) {
        std::ofstream outfile(file);
        outfile.setf(std::ios::scientific);
        outfile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        if (ndims == 1) {
            // Vector
            const size_t imax = reinterpret_cast<size_t*>(dims)[0];
            outfile << "%%MatrixMarket vector coordinate real general \n";
            for (size_t i = 0; i < ndims; ++i)
                outfile << reinterpret_cast<size_t*>(dims)[i] << " ";
            outfile << count_non_zero(*this) << std::endl;
            for (size_t i = 0; i < imax; ++i) {
                // NIST MatrixMarket version 3.0
                // Indexes are one-based
                const size_t idx1 = i + 1;
                if (this->vectorGet(i) != 0)
                    outfile << idx1 << "\t" << this->vectorGet(i) << std::endl;
            }
        } else if (ndims == 2) {
            // Matrix
            const size_t imax = reinterpret_cast<size_t*>(dims)[0];
            const size_t jmax = reinterpret_cast<size_t*>(dims)[1];
            outfile << "%%MatrixMarket matrix coordinate real general \n";
            for (size_t i = 0; i < ndims; ++i)
                outfile << reinterpret_cast<size_t*>(dims)[i] << " ";
            outfile << count_non_zero(*this) << std::endl;
            for (size_t i = 0; i < imax; ++i) {
                for (size_t j = 0; j < jmax; ++j) {
                    // NIST MatrixMarket version 3.0
                    // Indexes are one-based
                    const size_t idx1 = i + 1;
                    const size_t idx2 = j + 1;
                    if (this->matrixGet(i, j) != 0)
                        outfile << idx1 << "\t" << idx2 << "\t" << this->matrixGet(i, j) << std::endl;
                }
            }
        } else if (ndims == 3) {
            // T3D
            const size_t imax = reinterpret_cast<size_t*>(dims)[0];
            const size_t jmax = reinterpret_cast<size_t*>(dims)[1];
            const size_t kmax = reinterpret_cast<size_t*>(dims)[2];
            outfile << "%%MatrixMarket t3d coordinate real general \n";
            for (size_t i = 0; i < ndims; ++i)
                outfile << reinterpret_cast<size_t*>(dims)[i] << " ";
            outfile << count_non_zero(*this) << std::endl;
            for (size_t i = 0; i < imax; ++i) {
                for (size_t j = 0; j < jmax; ++j) {
                    for (size_t k = 0; k < kmax; ++k) {
                        const size_t idx1 = i + 1;
                        const size_t idx2 = j + 1;
                        const size_t idx3 = k + 1;
                        if (this->tensor3DGet(i, j, k) != 0)
                            outfile << idx1 << "\t" << idx2 << "\t" << idx3 << "\t" << this->tensor3DGet(i, j, k) << std::endl;
                    }
                }
            }
        } else if (ndims == 4) {
            const size_t imax = reinterpret_cast<size_t*>(dims)[0];
            const size_t jmax = reinterpret_cast<size_t*>(dims)[1];
            const size_t kmax = reinterpret_cast<size_t*>(dims)[2];
            const size_t lmax = reinterpret_cast<size_t*>(dims)[3];
            outfile << "%%MatrixMarket t4d coordinate real general \n";
            for (size_t i = 0; i < ndims; ++i)
                outfile << reinterpret_cast<size_t*>(dims)[i] << " ";
            outfile << count_non_zero(*this) << std::endl;
            for (size_t i = 0; i < imax; ++i) {
                for (size_t j = 0; j < jmax; ++j) {
                    for (size_t k = 0; k < kmax; ++k) {
                        for (size_t l = 0; l < lmax; ++l) {
                            const size_t idx1 = i + 1;
                            const size_t idx2 = j + 1;
                            const size_t idx3 = k + 1;
                            const size_t idx4 = l + 1;
                            if (this->tensor4DGet(i, j, k, l) != 0)
                                outfile << idx1 << "\t" << idx2 << "\t" << idx3 << "\t" << idx4 << "\t" << this->tensor4DGet(i, j, k, l) << std::endl;
                        }
                    }
                }
            }
        } else {
            HANDLE_ERROR("Unsupported ndims", 100);
        }
        outfile.close();
    }

    template<>
    template<>
    void LACIS<DENSE_LAPACK, double>::load_from_file<MATRIXMARKET>(const std::string& file) {
        std::ifstream infile(file);
        if (!infile.is_open()) {
            HANDLE_ERROR("Failed to open file: " + file, 101);
        }

        std::string line;
        std::getline(infile, line);
        if (line.find("%%MatrixMarket") != 0) {
            HANDLE_ERROR("Invalid MatrixMarket header", 102);
        }

        // Skip comments
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '%') continue;
            break;
        }

        std::istringstream dims_stream(line);
        std::vector<size_t> dims_vec;
        size_t val;
        while (dims_stream >> val) dims_vec.push_back(val);

        ndims = dims_vec.size() - 1;
        dims = new size_t[ndims];
        std::copy(dims_vec.begin(), dims_vec.end(), reinterpret_cast<size_t*>(dims));

        size_t nnz;
        if (!(infile >> nnz)) {
            HANDLE_ERROR("Missing number of non-zero entries", 103);
        }

        // Allocate memory
        if (ndims == 1)
            this->resizeToVector(reinterpret_cast<size_t*>(dims)[0]);
        else if (ndims == 2)
            this->resizeToMatrix(reinterpret_cast<size_t*>(dims)[0] , reinterpret_cast<size_t*>(dims)[1]);
        else if (ndims == 3)
            this->resizeToTensor3D(reinterpret_cast<size_t*>(dims)[0], reinterpret_cast<size_t*>(dims)[1], reinterpret_cast<size_t*>(dims)[2]);
        else if (ndims == 4)
            this->resizeToTensor4D(reinterpret_cast<size_t*>(dims)[0], reinterpret_cast<size_t*>(dims)[1], reinterpret_cast<size_t*>(dims)[2], reinterpret_cast<size_t*>(dims)[3]);
        else
            HANDLE_ERROR("Unsupported ndims", 102);

        for (size_t i = 0; i < nnz; ++i) {
            if (!std::getline(infile, line)) {
                HANDLE_ERROR("Unexpected EOF when reading values", 104);
            }
            std::istringstream iss(line);
            size_t i1, i2, i3, i4;
            double value;

            switch (ndims) {
                case 1:
                    if (!(iss >> i1 >> value)) HANDLE_ERROR("Invalid vector entry", 105);
                    this->vectorSet(i1 - 1, value);
                    break;

                case 2:
                    if (!(iss >> i1 >> i2 >> value)) HANDLE_ERROR("Invalid matrix entry", 107);
                    this->matrixSet(i1 - 1, i2 - 1, value);
                    break;

                case 3:
                    if (!(iss >> i1 >> i2 >> i3 >> value)) HANDLE_ERROR("Invalid 3D tensor entry", 109);
                    this->tensor3DSet(i1 - 1, i2 - 1, i3 - 1, value);
                    break;

                case 4:
                    if (!(iss >> i1 >> i2 >> i3 >> i4 >> value)) HANDLE_ERROR("Invalid 4D tensor entry", 111);
                    this->tensor4DSet(i1 - 1, i2 - 1, i3 - 1, i4 - 1, value);
                    break;

                default:
                    HANDLE_ERROR("Unsupported ndims", 113);
            }
        }

        infile.close();
    }


}
