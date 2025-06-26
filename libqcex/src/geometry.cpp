#include "qcex.hpp"
#include "qcex_utils.hpp"

#include <fstream>
#include <algorithm>
#include <sstream>

void _load_molecule_internal_geom (std::string fn, double** space) {
    LOG (INFO, "Loading molecular geometry from file "+fn);
    std::ifstream stream (fn);
    int molecule_size = 0;
    std::string top;
    std::getline(stream, top);

    std::ifstream copy (fn);

    if (top != "" || !is_number(top)) {
        std::cout << top << "\n";
        molecule_size = std::stoi(top);
    } else {
        HANDLE_ERROR ("Invalid header in file "+fn, 101);
    }

    long membytes = sizeof(double) * (4 * molecule_size + 1);

    LOG (DEV_INFO, "(malloc) Allocating geometry workspace (bytes): "+std::to_string(membytes));

    *space = reinterpret_cast<double*> (malloc(membytes));

    (*space)[0] = molecule_size; 
    int offset = 1;

    for (int i = 1; i <= molecule_size;) {
        std::string line;
        if (stream.eof() && (i - molecule_size) != 0) {
            LOG (DEV_INFO, "Expected more atoms: "+std::to_string(i - molecule_size));
            HANDLE_ERROR ("Reached EOF of geometry file", 103);
        } 
        std::getline(stream, line);
        if(line.empty()) continue;
        std::string symbol;
        double x, y, z;
        std::istringstream iss(line);
        if (!(stream >> symbol >> x >> y >> z)) {
            LOG (INFO, "Invalid parse: "+line);
            HANDLE_ERROR ("Could not parse geometry file", 104);
        }
        int j = i - 1;
        symbol = normalize_symbol(symbol);
        (*space)[4 * j + 0 + offset] = get_atomic_number (symbol);
        (*space)[4 * j + 1 + offset] = x;
        (*space)[4 * j + 2 + offset] = y;
        (*space)[4 * j + 3 + offset] = z;
        i++;
    }

    LOG (INFO, "Geometry read and loaded into memory");
    stream.close();
}

namespace newscf::qcex {

    void load_molecule (std::string filename, double** molecule) {
        _load_molecule_internal_geom (filename, molecule);
    }

    void destroy_molecule (double* molecule) {
        LOG (DEV_INFO, "(free) Removing geometry workspace");
        free(molecule);
    }

}