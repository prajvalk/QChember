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

#include <newscf/cis.hpp>
#include <newscf/newscf.hpp>
#include <newscf/newscf_utils.hpp>

#include <fstream>
#include <sstream>
#include <set>

namespace newscf {

    int load_basis (const std::string& filename, double** basis) {
		using namespace newscf::cis;
        std::ifstream infile(filename);
	    if (!infile.is_open()) {
	        HANDLE_ERROR("Failed to open basis file " + filename, 105);
	    }

    	// ignore line 1
    	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	    // === FIRST PASS TO CALCULATE MEMORY ===
	    int ndoubles = 0;
	    bool reading_atom = false;
	    std::string line;

	    while (std::getline(infile, line)) {
	    	if (utils::starts_with(line, "!")) continue;
	        if (utils::starts_with(line, "****")) {
	            reading_atom = true;
	            ndoubles += 2; // Atom Z, NumShells
	            continue;
	        }

	        if (reading_atom) {
	            std::istringstream iss(line);
	            std::string shell;
	            int nprim;
	            double scale;
	            if ((iss >> shell >> nprim >> scale) && std::isalpha(shell[0])) {
	                if (nprim <= 0) {
	                    HANDLE_ERROR("Invalid number of primitives <= 0 in basis: " + std::to_string(nprim), 99002);
	                }
	                ndoubles += 2 + 2 * nprim;  // L, nprim, (exp, coef) * nprim
	            }
	        }
	    }

	    infile.close();

	    size_t membytes = ndoubles * sizeof(double);
    	void* memspace = MALLOC(membytes);
	    *basis = reinterpret_cast<double*>(memspace);
	    if (!*basis) {
	        HANDLE_ERROR("Memory allocation failed for basis", 107);
	    }

	    // === SECOND PASS TO LOAD DATA ===
    	infile.open(filename);
    	// ignore line 1
    	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    if (!infile.is_open()) {
	        HANDLE_ERROR("Failed to reopen basis file", 106);
	    }

	    size_t i = 0;
	    int current_Z = -1;
	    int shell_count = 0;
	    size_t shell_count_pos = 0;

	    while (std::getline(infile, line)) {
	    	if (utils::starts_with(line, "!")) continue;
	        if (utils::starts_with(line, "****")) {
	            if (current_Z != -1) {
	                (*basis)[shell_count_pos] = static_cast<double>(shell_count);
	            }
	            shell_count = 0;
	            current_Z = -1;
	            continue;
	        }

	        std::istringstream iss(line);
	        std::string token;
	        if (current_Z == -1 && (iss >> token)) {
	            current_Z = utils::get_atomic_number(token);
	            if (current_Z <= 0) {
	                HANDLE_ERROR("Invalid atomic number parsed: " + std::to_string(current_Z), 99003);
	            }
	            (*basis)[i++] = static_cast<double>(current_Z);
	            shell_count_pos = i++;
	        } else if (std::isalpha(line[0])) {
	            std::string shell_type;
	            int nprim;
	            double scale;
	            iss.clear();
	            iss.str(line);
	            iss >> shell_type >> nprim >> scale;

	            if (shell_type.length() >= 2) {
	                HANDLE_ERROR("Compound shells not supported.", 99001);
	            }
	            int L = utils::get_angular_momentum_from_symbol(shell_type[0]);
	            if (L < 0 || L > 4) {  // Only support up to G shell (L=4)
	                HANDLE_ERROR("Unsupported angular momentum L=" + std::to_string(L), 99004);
	            }
	            if (nprim <= 0) {
	                HANDLE_ERROR("Invalid number of primitives <= 0: " + std::to_string(nprim), 99005);
	            }

	            (*basis)[i++] = static_cast<double>(L);
	            (*basis)[i++] = static_cast<double>(nprim);
	            shell_count++;

	            for (int j = 0; j < nprim; ++j) {
	                if (!std::getline(infile, line)) {
	                    HANDLE_ERROR("Unexpected EOF reading primitives", 99006);
	                }
	                std::istringstream datastream(line);
	                std::string exp_str, coef_str;
	                if (!(datastream >> exp_str >> coef_str)) {
	                    HANDLE_ERROR("Failed to parse exponent/coefficient line: " + line, 99007);
	                }

	                double exp_val = utils::parse_fortran_double(exp_str);
	                double coef_val = utils::parse_fortran_double(coef_str);

	                if (exp_val <= 0.0) {
	                    HANDLE_ERROR("Invalid exponent value <= 0: " + std::to_string(exp_val), 99008);
	                }
	                if (coef_val == 0.0) {
	                    LOG(WARNING, "Coefficient is zero for exponent: " + std::to_string(exp_val));
	                }

	                (*basis)[i++] = exp_val;
	                (*basis)[i++] = coef_val;
	            }
	        }
	    }

	    if (current_Z != -1) {
	        (*basis)[shell_count_pos] = static_cast<double>(shell_count);
	    }

	    infile.close();

	    return ndoubles;
    }

	int build_uacache(const double* geom, int** cache, int natoms) {
    	DEVLOG(DEV_INFO, "Building unique atom cache");

    	if (natoms <= 0 || geom == nullptr || cache == nullptr)
    		return 0;

    	*cache = reinterpret_cast<int*>(malloc(sizeof(int) * natoms));
    	if (!*cache) HANDLE_ERROR("Memory allocation failed", 9001);

		for (int i = 0; i < natoms; ++i)
			(*cache)[i] = static_cast<int>(geom[4 * i]);

    	return natoms;
    }

	int load_basis(const std::string& filename, double** basis,
                        const double* geom, int n) {
	    using namespace newscf::cis;

	    std::ifstream infile(filename);
	    if (!infile.is_open()) {
	        HANDLE_ERROR("Failed to open basis file " + filename, 105);
	    }

    	int* atom_Z_list;
		int num_atoms = build_uacache(geom, &atom_Z_list, n);

	    // === Collect required atomic numbers ===
	    std::set<int> required_Z;
	    for (int i = 0; i < num_atoms; ++i) {
	        if (atom_Z_list[i] <= 0) {
	            HANDLE_ERROR("Invalid atomic number in molecule: " + std::to_string(atom_Z_list[i]), 99009);
	        }
	        required_Z.insert(atom_Z_list[i]);
	    }

	    // === FIRST PASS: Count how much memory is needed ===
    	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    int ndoubles = 0;
	    std::string line;
	    bool reading_atom = false;
	    int current_Z = -1;

	    while (std::getline(infile, line)) {
	    	if (utils::starts_with(line, "!")) continue;
	        if (utils::starts_with(line, "****")) {
	            reading_atom = true;
	            current_Z = -1;
	            continue;
	        }

	        std::istringstream iss(line);
	        std::string token;
	        if (reading_atom && current_Z == -1 && (iss >> token)) {
	            current_Z = utils::get_atomic_number(token);
	            if (required_Z.find(current_Z) != required_Z.end()) {
	                ndoubles += 2; // Z, shell count placeholder
	            } else {
	                reading_atom = false;
	            }
	        } else if (reading_atom && std::isalpha(line[0])) {
	            std::string shell;
	            int nprim;
	            double scale;
	            iss.clear(); iss.str(line);
	            iss >> shell >> nprim >> scale;

	            if (nprim <= 0)
	                HANDLE_ERROR("Invalid nprim <= 0", 99010);

	            ndoubles += 2 + 2 * nprim;
	        }
	    }

	    infile.close();

	    // === Allocate memory ===
	    void* memspace = MALLOC(ndoubles * sizeof(double));
	    *basis = reinterpret_cast<double*>(memspace);
	    if (!*basis) {
	        HANDLE_ERROR("Memory allocation failed for filtered basis", 107);
	    }

	    // === SECOND PASS: Load actual data ===
	    infile.open(filename);
	    if (!infile.is_open()) {
	        HANDLE_ERROR("Failed to reopen basis file", 106);
	    }
    	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	    size_t i = 0;
	    current_Z = -1;
	    size_t shell_count_pos = 0;
	    int shell_count = 0;
	    bool load_this_atom = false;

	    while (std::getline(infile, line)) {
	    	if (utils::starts_with(line, "!")) continue;
	        if (utils::starts_with(line, "****")) {
	            if (load_this_atom) {
	                (*basis)[shell_count_pos] = static_cast<double>(shell_count);
	            }
	            shell_count = 0;
	            current_Z = -1;
	            load_this_atom = false;
	            continue;
	        }

	        std::istringstream iss(line);
	        std::string token;
	        if (current_Z == -1 && (iss >> token)) {
	            current_Z = utils::get_atomic_number(token);
	            load_this_atom = (required_Z.find(current_Z) != required_Z.end());
	            if (load_this_atom) {
	                (*basis)[i++] = static_cast<double>(current_Z);
	                shell_count_pos = i++;
	            }
	        } else if (load_this_atom && std::isalpha(line[0])) {
	            std::string shell;
	            int nprim;
	            double scale;
	            iss.clear(); iss.str(line);
	            iss >> shell >> nprim >> scale;

	            if (shell.length() >= 2)
	                HANDLE_ERROR("Compound shells not supported", 99001);

	            int L = utils::get_angular_momentum_from_symbol(shell[0]);
	            if (L < 0 || L > 4)
	                HANDLE_ERROR("Unsupported L = " + std::to_string(L), 99004);

	            (*basis)[i++] = static_cast<double>(L);
	            (*basis)[i++] = static_cast<double>(nprim);
	            shell_count++;

	            for (int j = 0; j < nprim; ++j) {
	                if (!std::getline(infile, line))
	                    HANDLE_ERROR("Unexpected EOF in primitives", 99006);
	                std::istringstream prims(line);
	                std::string exp_str, coef_str;
	                prims >> exp_str >> coef_str;

	                double exp = utils::parse_fortran_double(exp_str);
	                double coef = utils::parse_fortran_double(coef_str);

	                if (exp <= 0.0)
	                    HANDLE_ERROR("Invalid exponent <= 0: " + std::to_string(exp), 99008);

	                if (coef == 0.0)
	                    LOG(WARNING, "Zero coef for exponent " + std::to_string(exp));

	                (*basis)[i++] = exp;
	                (*basis)[i++] = coef;
	            }
	        }
	    }

	    if (load_this_atom) {
	        (*basis)[shell_count_pos] = static_cast<double>(shell_count);
	    }

	    infile.close();
	    return static_cast<int>(i);
	}


    void destroy_basis (double* basis) {
    	using namespace newscf::cis;
        FREE(basis);
    }

}