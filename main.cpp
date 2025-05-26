#include "fr_solver.h"
#include "mesh.h"

int main()
{
    Mesh mesh{};
    init_fws_mesh(mesh, 10, 10, 100.0, 100.0, 30, 20);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 0.0, 0.0, ONE_STD_ATM});
    solver.set_fws_bc(AIR_RHO, 3 * MACH, ONE_STD_ATM);
    for (int i = 0; i < 10; ++i)
    {
        solver.advance(0.01);
    }
    return 0;
}
