#include "geomopt.hpp"

#include <iostream>

#include "newscf/cis.hpp"

#include <cmath>
#include <cassert>

namespace newscf {

	template<>
	void optimzeGeometry<PSO>(const OptimizationParameters &params) {
		using namespace newscf::cis;

		std::cout << "PSO optimization parameters:" << std::endl;
		std::cout << "Initial Geometry: " << params.geometry_size << " atoms" << std::endl;

		if (params.geometry == nullptr) {
			HANDLE_ERROR("Unset Geometry", 80);
			return;
		}

		for (size_t i = 0; i < params.geometry_size; i++)
			std::cout << "\t" << params.geometry[i] << std::endl;

		std::cout << std::endl;

		std::cout << "Initial Search Space: " << params.searchspace_size << " parameters" << std::endl;

		if (params.searchspace_bounds == nullptr) {
			HANDLE_ERROR("Unset Searchspace", 81);
			return;
		}

		for (size_t i = 0; i < params.searchspace_size; i++)
			std::cout << "\t [" << params.searchspace_bounds[2*i] << "," << params.searchspace_bounds[2*i+1] << "]" << std::endl;

		std::cout << std::endl;

		std::cout << "Swarm Particles: " << params.pso_num_particles << std::endl;
		std::cout << "Swarm Inertial Weight: " << params.pso_cognitive_weight << std::endl;
		std::cout << "Swarm Social Weight: " << params.pso_social_weight << std::endl;

		typedef double psonum_t;

		// --- 1) compute element counts (all element counts, not bytes)
		const size_t sp = params.searchspace_size;
		const size_t n  = params.pso_num_particles;

		// sizes (in number of psonum_t elements)
		const size_t n_vbound_elems  = 2 * sp;       // lower/upper per dim for velocity bounds
		const size_t n_pos_elems     = n * sp;       // positions: particle-major, contiguous per particle
		const size_t n_vel_elems     = n * sp;       // velocities: same layout as positions
		const size_t n_pbestvec_elems= n * sp;       // personal best vectors
		const size_t n_gbestvec_elems= sp;           // global best vector
		const size_t n_pbestval_elems= n;            // per-particle best values
		const size_t n_gbestval_elems= 1;

		// total elements
		const size_t total_elems = n_vbound_elems + n_pos_elems + n_vel_elems +
		                     n_pbestvec_elems + n_gbestvec_elems +
		                     n_pbestval_elems + n_gbestval_elems;

		// bytes
		const size_t workspace_size = total_elems * sizeof(psonum_t);

		// offsets (in elements)
		const size_t OFF_VBOUNDS   = 0;
		const size_t OFF_POSVEC    = OFF_VBOUNDS   + n_vbound_elems;
		const size_t OFF_VELVEC    = OFF_POSVEC    + n_pos_elems;
		const size_t OFF_PBESTV    = OFF_VELVEC    + n_vel_elems;
		const size_t OFF_GBESTV    = OFF_PBESTV    + n_pbestvec_elems;
		const size_t OFF_PBESTVAL  = OFF_GBESTV    + n_gbestvec_elems;
		const size_t OFF_GBESTVAL  = OFF_PBESTVAL  + n_pbestval_elems;

		// allocate
		void* workspace = MALLOC(workspace_size);
		psonum_t* base = reinterpret_cast<psonum_t*>(workspace);

		psonum_t* VBOUNDS   = base + OFF_VBOUNDS;   // length = 2*sp, layout [low0, high0, low1, high1, ...]
		psonum_t* PPOS_V    = base + OFF_POSVEC;    // length = n*sp, layout particle-major: p0_dim0, p0_dim1, ... p1_dim0, ...
		psonum_t* PVEL_V    = base + OFF_VELVEC;    // same layout
		psonum_t* PBEST_V   = base + OFF_PBESTV;    // same layout
		psonum_t* GBEST_V   = base + OFF_GBESTV;    // length = sp
		psonum_t* PBEST_VAL = base + OFF_PBESTVAL;  // length = n
		psonum_t* GBEST_VAL = base + OFF_GBESTVAL;  // single element

		// sanity
		assert((OFF_GBESTVAL + n_gbestval_elems) == total_elems);

		// --- 2) build symmetric velocity bounds per-dimension
		for (size_t j = 0; j < sp; ++j) {
		    const psonum_t low  = params.searchspace_bounds[2*j + 0];
		    const psonum_t high = params.searchspace_bounds[2*j + 1];
		    const psonum_t span = fabs(high - low);
		    const psonum_t sp_width = params.pso_velocity_scaling_factor * span;

		    VBOUNDS[2*j + 0] = -sp_width;   // vmin
		    VBOUNDS[2*j + 1] = +sp_width;   // vmax
		}

		// --- 3) random generators -- pass VBOUNDS now (they are final)
		CISRandomGenerator<float> random_generator(512, 0, 1);
		CISRandomGenerator2<psonum_t> uniform_velocity_generator (sp, VBOUNDS);
		CISRandomGenerator2<psonum_t> uniform_search_space_generator (sp, params.searchspace_bounds);

		// --- init
		*GBEST_VAL = std::numeric_limits<psonum_t>::infinity();
		for (size_t i = 0; i < n; ++i) {
		    psonum_t* pos = PPOS_V + i*sp;
		    psonum_t* vel = PVEL_V + i*sp;
		    psonum_t* pbestvec = PBEST_V + i*sp;

		    uniform_search_space_generator.next(pos);  // fills sp values per particle
		    uniform_velocity_generator.next(vel);      // fills sp velocities

		    PBEST_VAL[i] = params.objective(pos, sp);
		    MEMCPY(pbestvec, pos, sp * sizeof(psonum_t));

		    if (PBEST_VAL[i] < *GBEST_VAL) {
		        *GBEST_VAL = PBEST_VAL[i];
		        MEMCPY(GBEST_V, pos, sp * sizeof(psonum_t));
		    }
		}

		// --- 4) main velocity-update loop
		const double w = params.pso_max_inertia;
		const double c1 = params.pso_cognitive_weight;
		const double c2 = params.pso_social_weight;

		for (size_t iter = 0; iter <= params.pso_max_iterations; iter++) {
			for (size_t i = 0; i < n; ++i) {
				psonum_t* pos = PPOS_V + i*sp;
				psonum_t* vel = PVEL_V + i*sp;
				psonum_t* pbest = PBEST_V + i*sp;

				for (size_t j = 0; j < sp; ++j) {
					const double r1 = random_generator.next();
					const double r2 = random_generator.next();

					// velocity update with cognitive/social coefficients
					double dv = w * vel[j]
								+ c1 * r1 * (pbest[j] - pos[j])
								+ c2 * r2 * (GBEST_V[j] - pos[j]);

					// clamp velocity
					psonum_t vmin = VBOUNDS[2 * j + 0];
					psonum_t vmax = VBOUNDS[2 * j + 1];
					if (dv < vmin) dv = vmin;
					if (dv > vmax) dv = vmax;
					vel[j] = dv;

					// position update
					pos[j] += vel[j];

					// position bound handling: reflect into bounds (or clamp)
					psonum_t xmin = params.searchspace_bounds[2 * j + 0];
					psonum_t xmax = params.searchspace_bounds[2 * j + 1];
					if (pos[j] < xmin) {
						pos[j] = xmin + (xmin - pos[j]);    // reflection
						if (pos[j] > xmax) pos[j] = xmin;  // fallback clamp
						vel[j] = 0; // damp velocity after hitting boundary
					} else if (pos[j] > xmax) {
						pos[j] = xmax - (pos[j] - xmax);    // reflection
						if (pos[j] < xmin) pos[j] = xmax;
						vel[j] = 0;
					}

					// explore objective value
					const psonum_t objval = params.objective(pos, sp);

					if (objval < PBEST_VAL[i]) {
						PBEST_VAL[i] = objval;
						MEMCPY(pbest, pos, sp * sizeof(psonum_t));
					}

					if (objval < *GBEST_VAL) {
						*GBEST_VAL = objval;
						MEMCPY(GBEST_V, pos, sp * sizeof(psonum_t));
					}
				}
			}
			std::cout << "Iteration #" << (iter + 1) << "\t Swarm Best: " << *GBEST_VAL << "\n";
		}
	}


}