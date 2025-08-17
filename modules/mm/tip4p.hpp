//
// Created by chlorinepentoxide on 7/8/25.
//

#ifndef TIP4P_HPP
#define TIP4P_HPP

#include <cmath>
#include "mm.hpp"

#if not defined USE_TIP4P_ORIGINAL_MODEL or not defined USE_TIP4P_ICE_MODEL or not defined USE_TIP4P_2005_MODEL or not defined USE_TIP4P_EWALD_MODEL
#define USE_TIP4P_ORIGINAL_MODEL
#endif

namespace newscf::mm::tip4p {
	namespace lj {

#ifdef USE_TIP4P_ORIGINAL_MODEL
		const double oo_epsilon = 0.1550; // kcal/mol
		const double oo_sigma   = 3.1536; // ang
#endif

#ifdef USE_TIP4P_ICE_MODEL
		const double oo_epsilon = 0.21084; // kcal/mol
		const double oo_sigma   = 3.1668; // ang
#endif

#ifdef USE_TIP4P_2005_MODEL
		const double oo_epsilon = 0.1852; // kcal/mol
		const double oo_sigma   = 3.1589; // ang
#endif

#ifdef USE_TIP4P_EWALD_MODEL
		const double oo_epsilon = 0.16275; // kcal/mol
		const double oo_sigma   = 3.16435; // ang
#endif

	}

	namespace charges {
#ifdef USE_TIP4P_ORIGINAL_MODEL
		const double o_charge = -1.040; // e
		const double h_charge = 0.520; // e
#endif

#ifdef USE_TIP4P_ICE_MODEL
		const double o_charge = -1.1794; // e
		const double h_charge = 0.5897; // e
#endif

#ifdef USE_TIP4P_2005_MODEL
		const double o_charge = -1.1128; // e
		const double h_charge = 0.5564; // e
#endif

#ifdef USE_TIP4P_EWALD_MODEL
		const double o_charge = -1.04844; // e
		const double h_charge = 0.52422; // e
#endif
	}

	namespace structure {
		const double r0_oh_bond = 0.9572; // ang
		const double theta0_hoh_bond = 104.52 * (M_PI / 180); // rad

#ifdef USE_TIP4P_ORIGINAL_MODEL
		const double om_distance = 0.15; // ang
#endif

#ifdef USE_TIP4P_ICE_MODEL
		const double om_distance = 0.1577; // ang
#endif

#ifdef USE_TIP4P_2005_MODEL
		const double om_distance = 0.1546; // ang
#endif

#ifdef USE_TIP4P_EWALD_MODEL
		const double om_distance = 0.1250; // ang
#endif

		// Coordinates in Lab Frame
		inline const double* O  = new double[3] {0, 0, 0};
		inline const double* M  = new double[3] {om_distance, 0, 0};
		inline const double* H1 = new double[3] {r0_oh_bond * cos(theta0_hoh_bond / 2), 0, r0_oh_bond * sin(theta0_hoh_bond / 2)};
		inline const double* H2 = new double[3] {r0_oh_bond * cos(theta0_hoh_bond / 2), 0, -r0_oh_bond * sin(theta0_hoh_bond / 2)};

		constexpr size_t H1_OFF = 0;
		constexpr size_t O_OFF  = 3;
		constexpr size_t H2_OFF = 6;
		constexpr size_t M_OFF  = 9;

		// Writes a new water molecule into the result with the defined TIP4P reference model
		// Result must have atleast sizeof(double) * 3 * 4
		// Memory Alignment H1-O-H2-M
		inline void new_water(double* result) {
			result[H1_OFF + 0] = H1[0];
			result[H1_OFF + 1] = H1[1];
			result[H1_OFF + 2] = H1[2];

			result[O_OFF  + 0] = O[0];
			result[O_OFF  + 1] = O[1];
			result[O_OFF  + 2] = O[2];

			result[H2_OFF + 0] = H2[0];
			result[H2_OFF + 1] = H2[1];
			result[H2_OFF + 2] = H2[2];

			result[M_OFF  + 0] = M[0];
			result[M_OFF  + 1] = M[1];
			result[M_OFF  + 2] = M[2];
		}

		// Calculates Interaction Energy between two molecules of water
		inline double calculate_pair_interaction_energy(double* mol1, double* mol2) {
			double energy = 0.0;

			// convention (first index -> internal coordinate) (second index -> external coordinate)
			const double h11_h12_r = l2_distance(mol1 + H1_OFF, mol2 + H1_OFF);
			const double h21_h22_r = l2_distance(mol1 + H2_OFF, mol2 + H2_OFF);
			const double h11_h22_r = l2_distance(mol1 + H1_OFF, mol2 + H2_OFF);
			const double h21_h12_r = l2_distance(mol1 + H2_OFF, mol2 + H1_OFF);

			const double oo_r = l2_distance(mol1 + O_OFF, mol2 + O_OFF);
			const double mm_r = l2_distance(mol1 + M_OFF, mol2 + M_OFF);

			const double m1_h21_r = l2_distance(mol1 + M_OFF, mol2 + H1_OFF);
			const double m1_h22_r = l2_distance(mol1 + M_OFF, mol2 + H2_OFF);
			const double m2_h11_r = l2_distance(mol2 + M_OFF, mol1 + H1_OFF);
			const double m2_h12_r = l2_distance(mol2 + M_OFF, mol1 + H2_OFF);

			// LJ Interaction
			energy += 4 * lj::oo_epsilon * ( pow(lj::oo_sigma / oo_r, 12) - pow(lj::oo_sigma / oo_r, 6) );

			// Electrostatics
			energy += coulomb_energy_prefactor * (
				charges::o_charge * charges::o_charge / mm_r +
				charges::o_charge * charges::h_charge * (1/m1_h21_r + 1/m1_h22_r + 1/m2_h11_r + 1/m2_h12_r) +
				charges::h_charge * charges::h_charge * (1/h11_h12_r + 1/h21_h22_r + 1/h11_h22_r + 1/h21_h12_r)
			);

			return energy * 4.184; // convert kcal/mol to kJ/mol
		}
	}
}

#endif //TIP4P_HPP
