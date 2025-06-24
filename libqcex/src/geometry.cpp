#include "qcex.hpp"
#include "qcex_utils.hpp"

#include <fstream>
#include <algorithm>
#include <sstream>

void _load_molecule_internal_geom (std::string fn, double* space) {
    LOG (INFO, "Loading molecular geometry from file "+fn);
    std::ifstream stream (fn);
    int molecule_size = 0;
    std::string top;
    std::getline(stream, top);

    try {
        if (top != "" || !is_number(top)) {
            molecule_size = std::stoi(top);
        } else {
            HANDLE_ERROR ("Invalid header in file "+fn, 101);
        }
    } catch(const std::exception& e) {
        HANDLE_EXCEPTION (e, "Internal error in processing geometry header", 102);
    }

    long membytes = sizeof(double) * 4 * molecule_size;

    LOG (DEV_INFO, "(malloc) Allocating geometry workspace (bytes): "+membytes);

    space = new double[molecule_size * 4];

    for (int i = 1; i <= molecule_size;) {
        std::string line;
        if (stream.eof() && (i - molecule_size) != 0) {
            //LOG (DEV_INFO, "Expected more atoms: "+(i - molecule_size));
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
        space[4 * j + 0] = get_atomic_number (symbol);
        space[4 * j + 1] = x;
        space[4 * j + 2] = y;
        space[4 * j + 3] = z;
        i++;
    }

    LOG (INFO, "Geometry read and loaded into memory");
    stream.close();
}

namespace newscf::qcex {

    void load_molecule (std::string filename, double* molecule) {
        _load_molecule_internal_geom (filename, molecule);
    }

    void destroy_molecule (double* molecule) {
        LOG (DEV_INFO, "(delete) Removing geometry workspace");
        delete[] molecule;
    }

}