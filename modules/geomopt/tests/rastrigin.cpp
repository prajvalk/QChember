//
// Created by chlorinepentoxide on 7/8/25.
//

#include <iostream>
#include <cmath>
#include "../geomopt.hpp"

double rastrigin_nd (const double* params, const int dim) {
	double result = 10 * dim;
	for (int i = 0; i < dim; ++i)
		result += (pow(params[i], 2) - cos(2 * M_PI * params[i]));
	return result;
}

int main() {
	int dim = 40 * 6;

	double* mol = new double[dim];
	double* searchspace = new double[dim * 2];

	for (int i = 0; i < dim; ++i) {
		mol[i] = pow(-1, i) * 0;
		searchspace[2 * i + 0] = -5.12;
		searchspace[2 * i + 1] = +5.12;
	}

	using namespace newscf;

	parametric_function obj = [&] (const double* params, int nparams) -> double {
		return rastrigin_nd(params, nparams);
	};

	OptimizationParameters params;
	params.searchspace_bounds = searchspace;
	params.searchspace_size = dim;
	params.geometry = mol;
	params.geometry_size = dim;
	params.objective = obj;

	optimzeGeometry<PSO>(params);

	std::cout << rastrigin_nd(mol, dim) << std::endl;

	delete[] mol;
}