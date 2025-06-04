#include <gtest/gtest.h>

#include "shape_f.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(MAPPING, face)
{
    auto bottom1 = face_mapping(BOTTOM);
    auto bottom2 = face_mapping((TOP + 2) % 4);
    EXPECT_EQ(bottom1, bottom2);
}

TEST(Lagrange, shape_func)
{
    EXPECT_EQ(lagrange(0, -1), 1);
    EXPECT_EQ(lagrange(1, -1), 0);
    EXPECT_EQ(lagrange(2, -1), 0);
    EXPECT_EQ(lagrange(0, 0), 0);
    EXPECT_EQ(lagrange(1, 0), 1);
    EXPECT_EQ(lagrange(2, 0), 0);
    EXPECT_EQ(lagrange(0, 1), 0);
    EXPECT_EQ(lagrange(1, 1), 0);
    EXPECT_EQ(lagrange(2, 1), 1);
}
