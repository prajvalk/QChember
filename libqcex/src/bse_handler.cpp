#include "qcex.hpp"

#include "logging_api.hpp"

void _bse_handler_load_basis (std::string filename, double** basis) {
    // Handles Psi4 basis format only
    
    
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