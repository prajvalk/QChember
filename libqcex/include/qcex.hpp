/*
 * Main API Header
 * new_scf Project
 * MPL-2.0 License
 */

#ifndef NEWSCF_QCEX_HPP
#define NEWSCF_QCEX_HPP

#include <string>
#include "matrix.hpp"

namespace newscf::qcex {

    using newscf::libmatrix::Matrix;

    // Level 0 Functions -- Optimized Low-Level Routines
    
    // Loads geometry of molecule into memory; Only XYZ Files are supported
    int  load_molecule    (std::string filename, double** molecule);
    // Destroys all memory used to store geometry
    void destroy_molecule (double* molecule);

    // Loads a complete basis from file to memory; Only Psi4 format is supported
    int  load_basis       (std::string filename, double** basis);
    // Loads a partial basis; only basis functions relating to the molecule information
    int  load_basis       (std::string filename, double** basis, double* molecule);
    // Destroys all memory associated with the basis
    void destroy_basis    (double* basis);

    // Level 1 Functions -- Wrapper Code

    struct IntegralEngineHandle {
        int*    atm;       // Atomic Description
        int*    bas;       // Basis Description
        double* env;       // Environment Description
        int*    shell_to_index; // Internal Shell Indexing Map
        double* cache;     // libcint cache
        int     natm;      // Number of Atoms
        int     nbas;      // Number of Basis Shells
        int     nenv;
        int     nbf_tot;   // Number of Basis Functions (Total)
        int     nbf_max;   // Number of Basis Functions (Max in any shell) 
    };

    void  intialize_handle (double* molecule, int molecule_size, double* basis, int basis_size, IntegralEngineHandle* handle);
    void  destroy_handle   (IntegralEngineHandle* handle);

    void  calculate_overlap_matrix             (IntegralEngineHandle* handle, Matrix<double>** output);
    void  calculate_kinetic_energy_matrix      (IntegralEngineHandle* handle, Matrix<double>** output);
    void  calculate_nuclear_attraction_matrix  (IntegralEngineHandle* handle, Matrix<double>** output);

    // Level 2 Functions -- Quantum Chem Drivers

    // Options for Calculations
    enum ComputeType {
        DIRECT,      // Integral direct ERI Calculations
        AOT_INCORE,  // Ahead of Time, incore ERI Calculations
        AOT_DISK     // Ahead of Time, disk written Calculations
    };

    struct Options {
        ComputeType computeType;
    };
}

#endif