#include <gtest/gtest.h>

#include "corr_f.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(Corr_F, g2)
{
    EXPECT_DOUBLE_EQ(g2L(-1.0), 1.0);
    EXPECT_DOUBLE_EQ(g2L(1.0), 0.0);
    EXPECT_DOUBLE_EQ(g2R(-1.0), 0.0);
    EXPECT_DOUBLE_EQ(g2R(1.0), 1.0);
}
