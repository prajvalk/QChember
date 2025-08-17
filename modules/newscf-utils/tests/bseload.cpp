#include "newscf/newscf.hpp"
#include "newscf/testing_cis.hpp"

int main() {
    using namespace newscf::cis;
    using namespace newscf;
    using namespace newscf::testing_cis;

    TestHandle h;
    h.testsuite = "basis_load_test";

    init_tests();

    add_test("H/def2-SVP test (complete basis load)", h);

    std::string bas1 = NEWSCF_TEST_DIR;
    bas1 += "/basis_sets/def2-svp.1.gbs";
    double* bswork;
    load_basis(bas1, &bswork);

    // Atom Z and Shell Count
    TEST_ASSERT_EQ(bswork[0], 1.0, h); // Z value for H
    TEST_ASSERT_EQ(bswork[1], 3.0, h); // Number of shells

    // Shell 1: S, 3 primitives
    TEST_ASSERT_EQ(bswork[2], 0.0, h); // Angular momentum L = 0 (S)
    TEST_ASSERT_EQ(bswork[3], 3.0, h); // 3 primitives

    TEST_ASSERT_EQ(bswork[4], 13.010701, h);
    TEST_ASSERT_EQ(bswork[5], 0.019682158, h);
    TEST_ASSERT_EQ(bswork[6], 1.9622572, h);
    TEST_ASSERT_EQ(bswork[7], 0.13796524, h);
    TEST_ASSERT_EQ(bswork[8], 0.44453796, h);
    TEST_ASSERT_EQ(bswork[9], 0.47831935, h);

    // Shell 2: S, 1 primitive
    TEST_ASSERT_EQ(bswork[10], 0.0, h); // Angular momentum L = 0
    TEST_ASSERT_EQ(bswork[11], 1.0, h);

    TEST_ASSERT_EQ(bswork[12], 0.12194962, h);
    TEST_ASSERT_EQ(bswork[13], 1.0, h);

    // Shell 3: P, 1 primitive
    TEST_ASSERT_EQ(bswork[14], 1.0, h); // Angular momentum L = 1
    TEST_ASSERT_EQ(bswork[15], 1.0, h);

    TEST_ASSERT_EQ(bswork[16], 0.8, h);
    TEST_ASSERT_EQ(bswork[17], 1.0, h);

    destroy_basis(bswork);
    complete_test(h);

    add_test("Acetone/def2-TZVPP test (partial basis test)", h);
    std::string mol = NEWSCF_TEST_DIR;
    mol += "/geometries/testsuite_bronze/03_Methanol.xyz";
    double* methanol;
    const int nmol = load_molecule(mol, &methanol);

    std::string bas2 = NEWSCF_TEST_DIR;
    bas2 += "/basis_sets/def2-tzvpp.1.gbs";
    load_basis(bas2, &bswork, methanol, nmol);

        // === H ===
    TEST_ASSERT_EQ(bswork[0], 1.0, h);  // Z
    TEST_ASSERT_EQ(bswork[1], 6.0, h);  // NumShells

    TEST_ASSERT_EQ(bswork[2], 0.0, h); // L
    TEST_ASSERT_EQ(bswork[3], 3.0, h);
    TEST_ASSERT_EQ(bswork[4], 34.0613410, h);
    TEST_ASSERT_EQ(bswork[5], 0.0060251978, h);
    TEST_ASSERT_EQ(bswork[6], 5.1235746, h);
    TEST_ASSERT_EQ(bswork[7], 0.045021094, h);
    TEST_ASSERT_EQ(bswork[8], 1.1646626, h);
    TEST_ASSERT_EQ(bswork[9], 0.20189726, h);

    TEST_ASSERT_EQ(bswork[10], 0.0, h);
    TEST_ASSERT_EQ(bswork[11], 1.0, h);
    TEST_ASSERT_EQ(bswork[12], 0.32723041, h);
    TEST_ASSERT_EQ(bswork[13], 1.0, h);

    TEST_ASSERT_EQ(bswork[14], 0.0, h);
    TEST_ASSERT_EQ(bswork[15], 1.0, h);
    TEST_ASSERT_EQ(bswork[16], 0.10307241, h);
    TEST_ASSERT_EQ(bswork[17], 1.0, h);

    TEST_ASSERT_EQ(bswork[18], 1.0, h);
    TEST_ASSERT_EQ(bswork[19], 1.0, h);
    TEST_ASSERT_EQ(bswork[20], 1.407, h);
    TEST_ASSERT_EQ(bswork[21], 1.0, h);

    TEST_ASSERT_EQ(bswork[22], 1.0, h);
    TEST_ASSERT_EQ(bswork[23], 1.0, h);
    TEST_ASSERT_EQ(bswork[24], 0.388, h);
    TEST_ASSERT_EQ(bswork[25], 1.0, h);

    TEST_ASSERT_EQ(bswork[26], 2.0, h);
    TEST_ASSERT_EQ(bswork[27], 1.0, h);
    TEST_ASSERT_EQ(bswork[28], 1.057, h);
    TEST_ASSERT_EQ(bswork[29], 1.0, h);

    // === C ===
    TEST_ASSERT_EQ(bswork[30], 6.0, h);
    TEST_ASSERT_EQ(bswork[31], 11.0, h);

    TEST_ASSERT_EQ(bswork[32], 0.0, h);
    TEST_ASSERT_EQ(bswork[33], 6.0, h);
    TEST_ASSERT_EQ(bswork[34], 13575.3496820, h);
    TEST_ASSERT_EQ(bswork[35], 0.00022245814352, h);
    TEST_ASSERT_EQ(bswork[36], 2035.2333680, h);
    TEST_ASSERT_EQ(bswork[37], 0.0017232738252, h);
    TEST_ASSERT_EQ(bswork[38], 463.22562359, h);
    TEST_ASSERT_EQ(bswork[39], 0.0089255715314, h);
    TEST_ASSERT_EQ(bswork[40], 131.20019598, h);
    TEST_ASSERT_EQ(bswork[41], 0.035727984502, h);
    TEST_ASSERT_EQ(bswork[42], 42.853015891, h);
    TEST_ASSERT_EQ(bswork[43], 0.11076259931, h);
    TEST_ASSERT_EQ(bswork[44], 15.584185766, h);
    TEST_ASSERT_EQ(bswork[45], 0.24295627626, h);

    TEST_ASSERT_EQ(bswork[46], 0.0, h);
    TEST_ASSERT_EQ(bswork[47], 2.0, h);
    TEST_ASSERT_EQ(bswork[48], 6.2067138508, h);
    TEST_ASSERT_EQ(bswork[49], 0.41440263448, h);
    TEST_ASSERT_EQ(bswork[50], 2.5764896527, h);
    TEST_ASSERT_EQ(bswork[51], 0.23744968655, h);

    TEST_ASSERT_EQ(bswork[52], 0.0, h);
    TEST_ASSERT_EQ(bswork[53], 1.0, h);
    TEST_ASSERT_EQ(bswork[54], 0.57696339419, h);
    TEST_ASSERT_EQ(bswork[55], 1.0, h);

    TEST_ASSERT_EQ(bswork[56], 0.0, h);
    TEST_ASSERT_EQ(bswork[57], 1.0, h);
    TEST_ASSERT_EQ(bswork[58], 0.22972831358, h);
    TEST_ASSERT_EQ(bswork[59], 1.0, h);

    TEST_ASSERT_EQ(bswork[60], 0.0, h);
    TEST_ASSERT_EQ(bswork[61], 1.0, h);
    TEST_ASSERT_EQ(bswork[62], 0.95164440028e-01 , h);
    TEST_ASSERT_EQ(bswork[63], 1.0, h);

    TEST_ASSERT_EQ(bswork[64], 1.0, h);
    TEST_ASSERT_EQ(bswork[65], 4.0, h);
    TEST_ASSERT_EQ(bswork[66], 34.697232244, h);
    TEST_ASSERT_EQ(bswork[67], 0.53333657805e-02, h);
    TEST_ASSERT_EQ(bswork[68], 7.9582622826, h);
    TEST_ASSERT_EQ(bswork[69], 0.035864109092, h);
    TEST_ASSERT_EQ(bswork[70], 2.3780826883, h);
    TEST_ASSERT_EQ(bswork[71], 0.14215873329, h);
    TEST_ASSERT_EQ(bswork[72], 0.81433208183, h);
    TEST_ASSERT_EQ(bswork[73], 0.34270471845, h);

    TEST_ASSERT_EQ(bswork[74], 1.0, h);
    TEST_ASSERT_EQ(bswork[75], 1.0, h);
    TEST_ASSERT_EQ(bswork[76], 0.28887547253, h);
    TEST_ASSERT_EQ(bswork[77], 0.46445822433, h);

    TEST_ASSERT_EQ(bswork[78], 1.0, h);
    TEST_ASSERT_EQ(bswork[79], 1.0, h);
    TEST_ASSERT_EQ(bswork[80], 0.10056823671, h);
    TEST_ASSERT_EQ(bswork[81], 0.24955789874, h);

    TEST_ASSERT_EQ(bswork[82], 2.0, h);
    TEST_ASSERT_EQ(bswork[83], 1.0, h);
    TEST_ASSERT_EQ(bswork[84], 1.097, h);
    TEST_ASSERT_EQ(bswork[85], 1.0, h);

    TEST_ASSERT_EQ(bswork[86], 2.0, h);
    TEST_ASSERT_EQ(bswork[87], 1.0, h);
    TEST_ASSERT_EQ(bswork[88], 0.318, h);
    TEST_ASSERT_EQ(bswork[89], 1.0, h);

    TEST_ASSERT_EQ(bswork[90], 3.0, h);
    TEST_ASSERT_EQ(bswork[91], 1.0, h);
    TEST_ASSERT_EQ(bswork[92], 0.761, h);
    TEST_ASSERT_EQ(bswork[93], 1.0, h);

    // === O ===
    TEST_ASSERT_EQ(bswork[94], 8.0, h);  // Z
    TEST_ASSERT_EQ(bswork[95], 11.0, h); // NumShells

    TEST_ASSERT_EQ(bswork[96], 0.0, h);  // L = 0 (S)
    TEST_ASSERT_EQ(bswork[97], 6.0, h);  // 6 primitives
    TEST_ASSERT_EQ(bswork[98], 27032.3826310, h);
    TEST_ASSERT_EQ(bswork[99], 0.00021726302465, h);
    TEST_ASSERT_EQ(bswork[100], 4052.3871392, h);
    TEST_ASSERT_EQ(bswork[101], 0.0016838662199, h);
    TEST_ASSERT_EQ(bswork[102], 922.32722710, h);
    TEST_ASSERT_EQ(bswork[103], 0.0087395616265, h);
    TEST_ASSERT_EQ(bswork[104], 261.24070989, h);
    TEST_ASSERT_EQ(bswork[105], 0.035239968808, h);
    TEST_ASSERT_EQ(bswork[106], 85.354641351, h);
    TEST_ASSERT_EQ(bswork[107], 0.11153519115, h);
    TEST_ASSERT_EQ(bswork[108], 31.035035245, h);
    TEST_ASSERT_EQ(bswork[109], 0.25588953961, h);

    TEST_ASSERT_EQ(bswork[110], 0.0, h); // S
    TEST_ASSERT_EQ(bswork[111], 2.0, h);
    TEST_ASSERT_EQ(bswork[112], 12.260860728, h);
    TEST_ASSERT_EQ(bswork[113], 0.39768730901, h);
    TEST_ASSERT_EQ(bswork[114], 4.9987076005, h);
    TEST_ASSERT_EQ(bswork[115], 0.24627849430, h);

    TEST_ASSERT_EQ(bswork[116], 0.0, h); // S
    TEST_ASSERT_EQ(bswork[117], 1.0, h);
    TEST_ASSERT_EQ(bswork[118], 1.1703108158, h);
    TEST_ASSERT_EQ(bswork[119], 1.0, h);

    TEST_ASSERT_EQ(bswork[120], 0.0, h); // S
    TEST_ASSERT_EQ(bswork[121], 1.0, h);
    TEST_ASSERT_EQ(bswork[122], 0.46474740994, h);
    TEST_ASSERT_EQ(bswork[123], 1.0, h);

    TEST_ASSERT_EQ(bswork[124], 0.0, h); // S
    TEST_ASSERT_EQ(bswork[125], 1.0, h);
    TEST_ASSERT_EQ(bswork[126], 0.18504536357, h);
    TEST_ASSERT_EQ(bswork[127], 1.0, h);

    TEST_ASSERT_EQ(bswork[128], 1.0, h); // P
    TEST_ASSERT_EQ(bswork[129], 4.0, h);
    TEST_ASSERT_EQ(bswork[130], 63.274954801, h);
    TEST_ASSERT_EQ(bswork[131], 0.60685103418e-02, h);
    TEST_ASSERT_EQ(bswork[132], 14.627049379, h);
    TEST_ASSERT_EQ(bswork[133], 0.041912575824, h);
    TEST_ASSERT_EQ(bswork[134], 4.4501223456, h);
    TEST_ASSERT_EQ(bswork[135], 0.16153841088, h);
    TEST_ASSERT_EQ(bswork[136], 1.5275799647, h);
    TEST_ASSERT_EQ(bswork[137], 0.35706951311, h);

    TEST_ASSERT_EQ(bswork[138], 1.0, h); // P
    TEST_ASSERT_EQ(bswork[139], 1.0, h);
    TEST_ASSERT_EQ(bswork[140], 0.52935117943, h);
    TEST_ASSERT_EQ(bswork[141], 0.44794207502, h);

    TEST_ASSERT_EQ(bswork[142], 1.0, h); // P
    TEST_ASSERT_EQ(bswork[143], 1.0, h);
    TEST_ASSERT_EQ(bswork[144], 0.17478421270, h);
    TEST_ASSERT_EQ(bswork[145], 0.24446069663, h);

    TEST_ASSERT_EQ(bswork[146], 2.0, h); // D
    TEST_ASSERT_EQ(bswork[147], 1.0, h);
    TEST_ASSERT_EQ(bswork[148], 2.314, h);
    TEST_ASSERT_EQ(bswork[149], 1.0, h);

    TEST_ASSERT_EQ(bswork[150], 2.0, h); // D
    TEST_ASSERT_EQ(bswork[151], 1.0, h);
    TEST_ASSERT_EQ(bswork[152], 0.645, h);
    TEST_ASSERT_EQ(bswork[153], 1.0, h);

    TEST_ASSERT_EQ(bswork[154], 3.0, h); // F
    TEST_ASSERT_EQ(bswork[155], 1.0, h);
    TEST_ASSERT_EQ(bswork[156], 1.428, h);
    TEST_ASSERT_EQ(bswork[157], 1.0, h);


    destroy_basis(bswork);

    end_tests(h);

    return h.exitcode;
}