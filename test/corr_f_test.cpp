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
    EXPECT_DOUBLE_EQ(g2L(0.5), g2R(-0.5));
}

TEST(Corr_F, dg2)
{
    EXPECT_DOUBLE_EQ(dg2L(-1.0), -dg2R(1.0));
    EXPECT_DOUBLE_EQ(dg2L(0.0), -dg2R(0.0));
    EXPECT_DOUBLE_EQ(dg2R(1.0), -dg2L(-1.0));
}
