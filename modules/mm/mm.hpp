//
// Created by chlorinepentoxide on 7/8/25.
//

#ifndef MM_HPP
#define MM_HPP

#include <cmath>

namespace newscf::mm {

    constexpr double coulomb_energy_prefactor = 3.32063711e+2L; // kcal/mol * Ang/e^2

	inline double l2_distance(const double* mol1, const double* mol2) {
		return sqrt(pow(mol1[0] - mol2[0], 2) + pow(mol1[1] - mol2[1], 2) + pow(mol1[2] - mol2[2], 2));
	}

}

#endif //MM_HPP
