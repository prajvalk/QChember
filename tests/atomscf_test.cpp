#include "newscf/newscf.hpp"
#include "newscf/common.hpp"

int main() {
	using namespace newscf;

	TestHandle h;
	h.testsuite = "atomscf";
	init_tests();

	add_test("He/def2-SVP", h);

	double* molecule;
	int mol_size = load_molecule("data/Helium.xyz", &molecule);

	double* basis;
	int basis_size = load_basis("data/def2-SVP.dat", &basis);

	IntegralEngineHandle inth;
	initialize_handle (molecule, mol_size, basis, basis_size, &inth);

	SolverOptions options;
	HFResult result;
	options.HF_STRONG_CONV = true;
	options.VERBOSE = true;
	atomscf_rhf (options, &result, &inth);

	TEST_ASSERT_EQ_TOL(result.HF_ENER, -2.855160479346704, options.HF_CONV_ETOL, h);

	complete_test(h);

	destroy_handle(&inth);
	destroy_basis(basis);

	add_test("He/aug-cc-pVQZ", h);

	basis_size = load_basis("data/aug-cc-pVQZ.dat", &basis);

	initialize_handle(molecule, mol_size, basis, basis_size, &inth);

	atomscf_rhf (options, &result, &inth);

	TEST_ASSERT_EQ_TOL(result.HF_ENER, -2.86152199563245, options.HF_CONV_ETOL, h);

	complete_test(h);

	destroy_handle (&inth);
	destroy_molecule(molecule);
	destroy_basis(basis);

	add_test("Ne/def2-SVP", h);

	mol_size = load_molecule("data/Neon.xyz", &molecule);
	basis_size = load_basis("data/def2-SVP.dat", &basis);
	initialize_handle(molecule, mol_size, basis, basis_size, &inth);
	atomscf_rhf (options, &result, &inth);

	destroy_handle (&inth);
	destroy_molecule(molecule);
	destroy_basis(basis);

	complete_test(h);

	end_tests(h);
	return h.exitcode;
}