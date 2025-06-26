#ifndef NEWSCF_QCEX_HPP
#define NEWSCF_QCEX_HPP

#include <string>

namespace newscf::qcex {

    void load_molecule (std::string filename, double** molecule);

    void destroy_molecule (double* molecule);

    void load_basis (std::string filename, double** basis);

    void load_basis (std::string filename, double** basis, double* molecule);

    void destroy_basis (double* basis);

}

#endif