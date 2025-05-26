#ifndef FR_SOLVER_H
#define FR_SOLVER_H

#include "euler_eq.h"
#include "type_def.h"
#include "shape_f.h"

// -----------------------------------------
// 核心类：FR 2D Euler 求解器
// -----------------------------------------

class FREulerSolver
{
public:
    FREulerSolver(Mesh mesh, const Vec4 &init_P);
    void setGamma(double g);
    void set_fws_bc(double rho, double u, double p);
    void advance(double dt);

private:
    const Mesh mesh;
    std::vector<Q9> nodes;
    double currentTime = 0.0;
    double gamma = GAMMA;
    std::array<double, 3> bc = {AIR_RHO, 1.0, ONE_STD_ATM}; // rho, u, p

    [[nodiscard]] std::vector<Q9> advance(const std::vector<Q9>& nodes, double dt) const;
    [[nodiscard]] Q9 computeResidual(int cellId) const;
};

#endif //FR_SOLVER_H
