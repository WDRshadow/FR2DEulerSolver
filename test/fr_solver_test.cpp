#include <gtest/gtest.h>

#include "fr_solver.h"
#include "mesh.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(Solve, runtime)
{
    Mesh mesh{};
    init_fws_mesh(mesh, 10, 10, 100.0, 100.0, 0.0, 0.0);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 3 * MACH, 0.0, ONE_STD_ATM});
    solver.advance(0.0001);
    EXPECT_EQ(1, 1);
}

TEST(Solve, vtu)
{
    Mesh mesh{};
    init_fws_mesh2(mesh, 2, 2, 1000.0, 1000.0, 1, 1);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 3 * MACH, 0.0, ONE_STD_ATM});
    solver.toVTU("output.vtu");
    EXPECT_EQ(1, 1);
}

TEST(Solve, convergence)
{
    Mesh mesh{};
    init_fws_mesh2(mesh, 10, 10, 1000.0, 1000.0, 7, 2);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 3 * MACH, 0.0, ONE_STD_ATM});
    for (int i = 0; i < 1000; ++i)
    {
        solver.advance(0.0001);
    }
    auto nodes = solver.getNodes();
    EXPECT_TRUE(nodes[2][8][0] > 0.0);
}
