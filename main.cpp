#include <iostream>

#include "fr_solver.h"
#include "mesh.h"

int main()
{
    Mesh mesh{};
    init_fws_mesh2(mesh, 30, 10, 3.0, 1.0, 24, 2);
    FREulerSolver solver(std::move(mesh), {1.4, 3.0, 0.0, 1.0});
    solver.set_fws_bc(1.4, 3.0, 0.0, 1.0);
    std::string fileName = "output/output_" + std::to_string(0) + ".vtu";
    solver.toVTU(fileName);
    for (int i = 1; i <= 1000; ++i)
    {
        solver.advance(0.0001);
        if (i % 10 == 0)
        {
            std::cout << "Current time: " << solver.getCurrentTime() << std::endl;
            fileName = "output/output_" + std::to_string(i / 10) + ".vtu";
            solver.toVTU(fileName);
        }
    }
    return 0;
}
