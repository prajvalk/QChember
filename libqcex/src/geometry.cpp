#include "qcex.hpp"
#include "qcex_utils.hpp"

#include <fstream>
#include <algorithm>
#include <sstream>

int _load_molecule_internal_geom(const std::string& fn, double** space) {
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

    long membytes = sizeof(double) * (1 + 4 * molecule_size);
    LOG(DEV_INFO, "(malloc) Allocating geometry workspace (bytes): " + std::to_string(membytes));
    *space = reinterpret_cast<double*>(malloc(membytes));
    if (!*space) {
        HANDLE_ERROR("Memory allocation failed for geometry", 104);
    }

    int offset = 0;

    int atoms_read = 0;
    while (atoms_read < molecule_size && std::getline(stream, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string symbol;
        double x, y, z;
        if (!(iss >> symbol >> x >> y >> z)) {
            HANDLE_ERROR("Failed to parse geometry line: " + line, 105);
        }

        symbol = normalize_symbol(symbol);
        int Z = get_atomic_number(symbol);

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

namespace newscf::qcex {

    int load_molecule (std::string filename, double** molecule) {
        return _load_molecule_internal_geom (filename, molecule);
    }

    void destroy_molecule (double* molecule) {
        LOG (DEV_INFO, "(free) Removing geometry workspace");
        free(molecule);
    }

}