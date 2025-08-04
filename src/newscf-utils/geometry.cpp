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
#include "newscf/newscf.hpp"
#include "newscf/newscf_utils.hpp"

#include <fstream>
#include <sstream>

namespace newscf {

    int load_molecule (const std::string& fn, double** space) {
		using namespace newscf::cis;

        LOG(INFO, "Loading molecular geometry from file " + fn);

		std::ifstream stream(fn);
		if (!stream.is_open()) {
		    HANDLE_ERROR("Failed to open file " + fn, 100);
		}

		std::string line;
		int molecule_size = 0;

		// === Read atom count line ===
		if (!std::getline(stream, line)) {
		    HANDLE_ERROR("Empty file or failed to read atom count line in " + fn, 101);
		}
		try {
		    molecule_size = std::stoi(line);
		    if (molecule_size <= 0) {
		        HANDLE_ERROR("Number of atoms must be positive, got: " + line, 102);
		    }
		} catch (const std::exception& e) {
		    HANDLE_ERROR("Failed to parse number of atoms from first line: " + line, 102);
		}

		// === Skip the comment/units line ===
		if (!std::getline(stream, line)) {
		    HANDLE_ERROR("Unexpected EOF after atom count line in " + fn, 103);
		}

		size_t membytes = sizeof(double) * (1 + 4 * molecule_size);
    	void* memspace = MALLOC(membytes);
		*space = reinterpret_cast<double*>(memspace);
		if (!*space) {
		    HANDLE_ERROR("Memory allocation failed for geometry", 104);
		}

		int atoms_read = 0;
		while (atoms_read < molecule_size && std::getline(stream, line)) {
			int offset = 0;
			if (line.empty()) continue;

		    std::istringstream iss(line);
		    std::string symbol;
		    double x, y, z;
		    if (!(iss >> symbol >> x >> y >> z)) {
		        HANDLE_ERROR("Failed to parse geometry line: " + line, 105);
		    }

		    symbol = utils::normalize_symbol(symbol);
		    int Z = utils::get_atomic_number(symbol);

		    size_t idx = offset + atoms_read * 4;
		    (*space)[idx + 0] = static_cast<double>(Z);
		    (*space)[idx + 1] = x;
		    (*space)[idx + 2] = y;
		    (*space)[idx + 3] = z;

		    atoms_read++;
		}

		if (atoms_read != molecule_size) {
		    HANDLE_ERROR("Expected " + std::to_string(molecule_size) +
		                 " atoms, but read " + std::to_string(atoms_read), 106);
		}

		LOG(INFO, "Successfully loaded " + std::to_string(molecule_size) + " atoms from geometry.");

		return molecule_size;
    }

    void destroy_molecule (double* molecule) {
        DEVLOG (cis::DEV_INFO, "(free) Removing geometry workspace");
        free(molecule);
    }

}