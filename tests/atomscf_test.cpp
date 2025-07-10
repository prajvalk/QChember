#include "newscf/newscf.hpp"
#include "newscf/common.hpp"

int main() {
	using namespace newscf;

	TestHandle h;
	h.testsuite = "atomscf";
	init_tests();

	double* molecule;
	int mol_size = load_molecule("data/Helium.xyz", &molecule);

	double* basis;
	int basis_size = load_basis("data/def2-SVP.dat", &basis);

	IntegralEngineHandle inth;
	initialize_handle (molecule, mol_size, basis, basis_size, &inth);

	SolverOptions options;
	HFResult result;
	options.HF_STRONG_CONV = true;
	atomscf_rhf (options, &result, &inth);

	TEST_ASSERT_EQ_TOL(result.HF_ENER, -2.855160479346704, options.HF_CONV_ETOL, h);

	complete_test(h);

	destroy_handle (&inth);
	destroy_molecule(molecule);
	destroy_basis(basis);

	end_tests(h);
	return h.exitcode;
}