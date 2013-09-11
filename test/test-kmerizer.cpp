#include "test.h"
#include "kmerizer/kmerizer.h"

Kmerizer * initKmerizer() {
    return new Kmerizer(1, 1, "/tmp", CANONICAL);
}

namespace {
    
class KmerizerTest : public ::testing::Test {
};

TEST(KmerizerTest, DoesInstantiation) {
    EXPECT_TRUE(initKmerizer() != NULL);
}

} /* namespace */

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
