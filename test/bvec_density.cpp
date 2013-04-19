#include "bvec/bvec.h"
#include "test.h"
#include <algorithm>

#define TEST_VEC_LENGTH 4
const uint32_t TEST_CASES[][TEST_VEC_LENGTH] = {
    { 161143029,  1089249216, 1093915745, 1157787156 },
    { 3085261341, 3247604874, 3456053010, 3773962658 },
    { 13431029,   377786868,  1463271703, 1501810865 },
    { 1819233732, 2235579467, 2625766209, 4061912144 }
};

using namespace std;

int run_test_cases() {
    int num_cases = sizeof(TEST_CASES) / (TEST_VEC_LENGTH * sizeof(uint32_t));
    int retValue = 0;
    for (int i = 0; i < num_cases; i++) {
        vector<uint32_t> v;
        for (int j = 0; j < TEST_VEC_LENGTH; j++)
            v.push_back(TEST_CASES[i][j]);
        bvec *bv = new bvec(v);
        bv->compress();
        bv->decompress();
        vector<uint32_t> uncomp = bv->get_words();
        bool result = v == uncomp;
        printf(""TERM_BOLD"Test %d"TERM_RESET": %s\n", i, result
            ? ""TERM_GREEN"Pass"TERM_RESET"" : ""TERM_RED"Fail"TERM_RESET"");
        if (!result) {
            vector<uint32_t> va;
            printf("%-20s ", "Original");
            for (int j = 0; j < v.size(); j++) printf("%32u ", v[j]);
            cout << endl;
            xor_vectors(v, uncomp, va);
            debug_binary("Original (bin)", v);
            debug_binary("Uncompressed", uncomp);
            debug_binary("XOR", va);
            retValue++;
        }
    }
    return retValue;
}

int main(int argc, char *argv[]) {
    return run_test_cases();
}