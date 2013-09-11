#include "test.h"
#include "bvec/bvec.h"
#include <algorithm>

#define TEST_VEC_LENGTH 4
const uint32_t TEST_CASES[][TEST_VEC_LENGTH] = {
    { 161143029,  1089249216, 1093915745, 1157787156 },
    { 3085261341, 3247604874, 3456053010, 3773962658 },
    { 13431029,   377786868,  1463271703, 1501810865 },
    { 1819233732, 2235579467, 2625766209, 4061912144 }
};

BitVector * random_bvec(int n) {
    // simulate random bitvectors with up to 1 billion points
    // 16mers - up to 2^32 - 1
    // 32mers - up to 2^64 - 1
    srand(time(NULL));

    vector<uint32_t> rand32;
    for(int i=0;i<n;i++) {
        rand32.push_back(rand() | ((rand() % 4) << 30));
    }
    sort(rand32.begin(),rand32.end());
    printf("max value selected: %u\n",rand32.back());
    vector<uint32_t>::iterator ii32 = rand32.begin();
    vector<uint32_t> uniq;
    uniq.push_back(*ii32);
    ++ii32;
    while (ii32 != rand32.end()) {
        if (*ii32 != uniq.back())
            uniq.push_back(*ii32);
        ++ii32;
    }
    printf("generated %zi distinct random uint32_t\n",uniq.size());
    return new BitVector(uniq);
}

namespace {
    
class BitVectorTest : public ::testing::Test {
};

TEST(BitVectorTest, CanDecompressCompressedObject) {
    int num_cases = sizeof(TEST_CASES) /
        (TEST_VEC_LENGTH * sizeof(uint32_t));
    int retValue = 0;
    for (int i = 0; i < num_cases; i++) {
        vector<uint32_t> v;
        for (int j = 0; j < TEST_VEC_LENGTH; j++)
            v.push_back(TEST_CASES[i][j]);
        BitVector *bv = new BitVector(v);
        bv->compress();
        bv->decompress();
        vector<uint32_t> uncomp = bv->getWords();
        EXPECT_EQ(v, uncomp);
    }
}

TEST(BitVectorTest, CanDeserializeSerializedObject) {
    short n = 1000;

    BitVector * original = random_bvec(n);
    BitVector * deserialized = new BitVector();

    char filename[200];
    std::tmpnam(filename);

    BitVector::save(*original, filename);
    BitVector::restore(*deserialized, filename);

    EXPECT_TRUE(original->equals(*deserialized));
}

} /* namespace */

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

/*
    bool match = original->equals(*deserialized); // original == deserialized
    if (!match) {
        original->print();
        deserialized->print();
        vector<uint32_t> orig = original->getWords();
        vector<uint32_t> deser = deserialized->getWords();
        debug_binary("Original", orig);
        debug_binary("Deserialized", deser);
        vector<uint32_t> v;
        xor_vectors(orig, deser, v);
        debug_binary("XOR", v);
        return 1;
    }
    return 0;
}
*/