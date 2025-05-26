#include "fr_solver.h"
#include "mesh.h"

int main()
{
    Mesh mesh{};
    init_fws_mesh(mesh, 10, 10, 1.0, 1.0, 0.0, 0.0);
    FREulerSolver solver(mesh, {AIR_RHO, 1.0, 0.0, ONE_STD_ATM});
    solver.advance(0.01);
    return 0;
}
