#include <iostream>

#include "fr_solver.h"
#include "mesh.h"

int main()
{
    Mesh mesh{};
    init_fws_mesh2(mesh, 10, 10, 1000.0, 1000.0, 7, 2);
    FREulerSolver solver(std::move(mesh), {AIR_RHO, 3 * MACH, 0.0, ONE_STD_ATM});
    std::string fileName = "output/output_" + std::to_string(0) + ".vtu";
    solver.toVTU(fileName);
    for (int i = 1; i <= 4000; ++i)
    {
        solver.advance(0.001);
        std::cout << "Current time: " << solver.getCurrentTime() << std::endl;
        fileName = "output/output_" + std::to_string(i) + ".vtu";
        solver.toVTU(fileName);
    }
    return 0;
}
