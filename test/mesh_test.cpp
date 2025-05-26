#include <gtest/gtest.h>

#include "mesh.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(Mesh, test)
{
    Mesh mesh{};
    init_fws_mesh(mesh, 2, 2, 1.0, 1.0, 0.5, 0.5);
    EXPECT_EQ(mesh.vertices.size(), 9);
    EXPECT_DOUBLE_EQ(mesh.vertices[4].y, 0.75);
    EXPECT_EQ(mesh.elements.size(), 4);
    EXPECT_EQ(mesh.elements[3].faceIds[2], 11);
    EXPECT_EQ(mesh.faces.size(), 12);
    EXPECT_EQ(mesh.faces[9].leftCell, 1);
}
