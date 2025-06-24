#ifndef NEWSCF_QCEX_HPP
#define NEWSCF_QCEX_HPP

#include <string>

#include "logging_api.hpp"

namespace newscf::qcex {

    void load_molecule (std::string filename, double* molecule);

    void destroy_molecule (double* molecule);

}

#endif