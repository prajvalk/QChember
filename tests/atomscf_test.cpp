#include "newscf/newscf.hpp"
#include "newscf/common.hpp"

int main() {
	using namespace newscf;

	TestHandle h;
	h.testsuite = "atomscf";
	init_tests();

	double* molecule;
	int mol_size = load_molecule("data/Water.xyz", &molecule);

	double* basis;
	int basis_size = load_basis("data/def2-SVP.dat", &basis);

	IntegralEngineHandle inth;
	initialize_handle (molecule, mol_size, basis, basis_size, &inth);

	SolverOptions options;
	atomscf_rhf (options, &inth);

	destroy_handle (&inth);
	destroy_molecule(molecule);
	destroy_basis(basis);

	end_tests(h);
	return h.exitcode;
}