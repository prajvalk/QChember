
#ifndef GEOMOPT_HPP
#define GEOMOPT_HPP
#include <cmath>
#include <functional>

namespace newscf {

	typedef std::function<double (const double* params, size_t nparams)> parametric_function;

	inline void euler_transform (const double* position, const double* transform, double* result) {
		const double x        = position[0];
		const double y        = position[1];
		const double z        = position[2];
		const double of_x     = transform[0];
		const double of_y     = transform[1];
		const double of_z     = transform[2];
		const double ro_theta = transform[3];
		const double ro_phi   = transform[4];
		const double ro_psi   = transform[5];

		result[0] = of_x + x * (cos(ro_psi) * cos(ro_phi) - sin(ro_psi) * sin(ro_phi) * cos(ro_theta))
		                 + y * (cos(ro_psi) * sin(ro_phi) + sin(ro_psi) * cos(ro_phi) * cos(ro_theta))
						 + z * (sin(ro_psi) * sin(ro_theta));

		result[1] = of_y + x * (-sin(ro_psi) * cos(ro_phi) - cos(ro_psi) * sin(ro_phi) * cos(ro_theta))
						 + y * (-sin(ro_psi) * sin(ro_phi) + cos(ro_psi) * cos(ro_phi) * cos(ro_theta))
						 + z * (cos(ro_psi) * sin(ro_theta));

		result[2] = of_z + x * sin(ro_theta) * sin(ro_phi)
						 + y * cos(ro_theta) * cos(ro_phi)
						 + z * cos(ro_theta);
	}

	enum GeometryOptimizers {
		PSO
	};

	struct OptimizationParameters {
		double* geometry;
		size_t geometry_size;
		double* searchspace_bounds;
		size_t searchspace_size;
		parametric_function objective;

		size_t pso_num_particles = 5000;
		double pso_max_inertia = 0.9;
		double pso_min_inertia = 0.4;
		double pso_cognitive_weight = 0.2;
		double pso_social_weight = 0.8;
		size_t pso_max_iterations = 500;
		double pso_velocity_scaling_factor = 0.2;
	};

	template<GeometryOptimizers opt>
	void optimzeGeometry(const OptimizationParameters& params);

}

#endif //GEOMOPT_HPP
