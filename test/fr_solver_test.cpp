#include <gtest/gtest.h>

#include "fr_solver.h"
#include "mesh.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(Solve, initialize)
{
    Mesh mesh{};
    init_fws_mesh(mesh, 10, 10, 1.0, 1.0, 0.0, 0.0);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 1.0, 0.0, ONE_STD_ATM});
    solver.advance(0.01);
    EXPECT_EQ(1, 1);
}

TEST(Solve, runtime)
{
    Mesh mesh{};
    init_fws_mesh(mesh, 10, 10, 100.0, 100.0, 30, 20);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 0.0, 0.0, ONE_STD_ATM});
    solver.set_fws_bc(AIR_RHO, 3 * MACH, ONE_STD_ATM);
    for (int i = 0; i < 10; ++i)
    {
        solver.advance(0.01);
    }
    EXPECT_EQ(1, 1);
}