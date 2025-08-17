#include <iostream>
#include <cmath>
#include "../geomopt.hpp"

double lj (const double* coords, const int nparams) {
	const int size = nparams / 3;
	double energy = 0;
	for (size_t i = 0; i < size; ++i) {
		const double x1 = coords[3 * i + 0];
		const double y1 = coords[3 * i + 1];
		const double z1 = coords[3 * i + 2];
		for (size_t j = 0; j < i; ++j) {
			const double x2 = coords[3 * j + 0];
			const double y2 = coords[3 * j + 1];
			const double z2 = coords[3 * j + 2];
			const double r = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
			energy += pow(1.0/r, 12) - pow(1.0/r, 6);
		}
	}
	return energy * 4 * 1;
}

int main() {
	const double bbox = 4;
	const int size = 10;

	double* ssp = new double[size * 3 * 2];
	for (int i = 0; i < 3 * size; ++i) {
		ssp[2*i + 0] = 0;
		ssp[2*i + 1] = bbox;
	}

	using namespace newscf;
	OptimizationParameters params;
	params.searchspace_bounds = ssp;
	params.searchspace_size = size * 3;

	double* mol = new double[1];
	mol[0] = 0;
	params.geometry = mol;
	params.geometry_size = 1;
	params.objective = lj;

	optimzeGeometry<PSO>(params);

	delete[] ssp;
	delete[] mol;
	return 0;
}