#include <fstream>
#include <sstream>

#include "qcex.hpp"
#include "qcex_utils.hpp"
#include "logging_api.hpp"
int _bse_handler_load_basis(std::string filename, double** basis) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        HANDLE_ERROR("Failed to open basis file " + filename, 105);
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
                if (nprim <= 0) {
                    HANDLE_ERROR("Invalid number of primitives <= 0 in basis: " + std::to_string(nprim), 99002);
                }
                ndoubles += 2 + 2 * nprim;  // L, nprim, (exp, coef) * nprim
            }
        }
    }

    infile.close();

    ndoubles += 1;  // For sentinel 0.0
    long membytes = ndoubles * sizeof(double);
    *basis = reinterpret_cast<double*>(malloc(membytes));
    if (!*basis) {
        HANDLE_ERROR("Memory allocation failed for basis", 107);
    }

    LOG(DEV_INFO, "(malloc) Allocating basis space (bytes): " + std::to_string(membytes));

    // === SECOND PASS TO LOAD DATA ===
    infile.open(filename);
    if (!infile.is_open()) {
        HANDLE_ERROR("Failed to reopen basis file", 106);
    }

    size_t i = 0;
    int current_Z = -1;
    int shell_count = 0;
    size_t shell_count_pos = 0;

    while (std::getline(infile, line)) {
        if (starts_with(line, "****")) {
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
            current_Z = get_atomic_number(token);
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
            int L = get_angular_momentum_from_symbol(shell_type[0]);
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

                double exp_val = parse_fortran_double(exp_str);
                double coef_val = parse_fortran_double(coef_str);

                if (exp_val <= 0.0) {
                    HANDLE_ERROR("Invalid exponent value <= 0: " + std::to_string(exp_val), 99008);
                }
                if (coef_val == 0.0) {
                    LOG(WARN, "Coefficient is zero for exponent: " + std::to_string(exp_val));
                }

                (*basis)[i++] = exp_val;
                (*basis)[i++] = coef_val;
            }
        }
    }

    if (current_Z != -1) {
        (*basis)[shell_count_pos] = static_cast<double>(shell_count);
    }

    // Write final sentinel
    (*basis)[i++] = 0.0;

    infile.close();

    return ndoubles;
}

namespace newscf::qcex {

    int load_basis (std::string filename, double** basis) {
        return _bse_handler_load_basis (filename, basis);
    }

    void destroy_basis (double* basis) {
        LOG (DEV_INFO, "(free) Removing basis workspace");
        free(basis);
    }

}