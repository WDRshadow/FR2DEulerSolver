#include <iostream>

#include "fr_solver.h"
#include "mesh.h"

int main()
{
    Mesh mesh{};
    init_inf_mesh(mesh, 10, 10, 10.0, 10.0);
    FREulerSolver solver(std::move(mesh), 0, 1, 5);
    std::string fileName = "output/output_" + std::to_string(0) + ".vtu";
    solver.toVTU(fileName);
    for (int i = 1; i <= 10000; ++i)
    {
        solver.advance(0.001);
        if (i % 10 == 0)
        {
            std::cout << "Current time: " << solver.getCurrentTime() << std::endl;
            fileName = "output/output_" + std::to_string(i / 10) + ".vtu";
            solver.toVTU(fileName);
        }
    }
    return 0;
}
