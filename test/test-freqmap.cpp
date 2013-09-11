#include "test.h"
#include "util/freqmap.h"

namespace {
    
class FrequencyMapTest : public ::testing::Test {
};

TEST(FrequencyMapTest, DoesInstantiation) {
    FrequencyMap *map = new FrequencyMap();
    EXPECT_TRUE(map != NULL);
}

TEST(FrequencyMapTest, CountsKmers) {
    FrequencyMap *map = new FrequencyMap();
    
    EXPECT_EQ(map->count("kmer1"), 0) << "Should be initialized to 0";
    EXPECT_EQ(map->count("kmer2"), 0) << "Should be initialized to 0";
    
    map->add("kmer1");
    map->add("kmer2");
    map->add("kmer1");
    
    EXPECT_EQ(map->count("kmer1"), 2) << "Should have a count of 2";
    EXPECT_EQ(map->count("kmer2"), 1) << "Should have a count of 1";
    EXPECT_EQ(map->count("kmer3"), 0) << "Should have a count of 0";
}

} /* namespace */

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
