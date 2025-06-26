#include <fstream>
#include <sstream>

#include "qcex.hpp"
#include "qcex_utils.hpp"
#include "logging_api.hpp"

void _bse_handler_load_basis (std::string filename, double** basis) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        HANDLE_ERROR ("Failed to open basis file "+filename, 105);
    }

    // === FIRST PASS TO CALCULATE MEMORY ===
    size_t ndoubles = 0;
    bool reading_atom = false;
    std::string line;

    while (std::getline(infile, line)) {
        if (starts_with(line, "****")) {
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
                ndoubles += 2 + 2 * nprim;  // L, nprim, (exp, coef) * nprim
            }
        }
    }

    infile.close();
    long membytes = ndoubles * sizeof(double);
    *basis = reinterpret_cast<double*>(malloc(membytes));

    // === SECOND PASS TO LOAD DATA ===
    infile.open(filename);
    if (!infile.is_open()) {
        HANDLE_ERROR ("Failed to reopen basis file", 106);
    }

    size_t i = 0;
    int current_Z = -1;
    int shell_count = 0;
    size_t shell_count_pos = 0;

    while (std::getline(infile, line)) {
        if (starts_with(line, "****")) {
            if (current_Z != -1) {
                *basis[shell_count_pos] = shell_count;
            }
            shell_count = 0;
            current_Z = -1;
            continue;
        }

        std::istringstream iss(line);
        std::string token;
        if (current_Z == -1 && (iss >> token)) {
            current_Z = get_atomic_number(token);
            *basis[i++] = current_Z;
            shell_count_pos = i++;
        } else if (std::isalpha(line[0])) {
            std::string shell_type;
            int nprim;
            double scale;
            iss.clear(); iss.str(line);
            iss >> shell_type >> nprim >> scale;
            if (shell_type.length() >= 2) {
                HANDLE_ERROR ("Compound shells not supported.", 99001);
            }
            int L = shell_type[0] - 'S';  // S=0, P=1, D=2, ...
            *basis[i++] = L;
            *basis[i++] = nprim;
            shell_count++;

            for (int j = 0; j < nprim; ++j) {
                std::getline(infile, line);
                std::istringstream datastream(line);
                double exp, coef;
                datastream >> exp >> coef;
                *basis[i++] = exp;
                *basis[i++] = coef;
            }
        }
    }

    if (current_Z != -1) {
        *basis[shell_count_pos] = shell_count;
    }

    infile.close();
}

namespace newscf::qcex {

    void load_basis (std::string filename, double** basis) {
        _bse_handler_load_basis (filename, basis);
    }

    void destroy_basis (double* basis) {
        LOG (DEV_INFO, "(free) Removing basis workspace");
        free(basis);
    }

}