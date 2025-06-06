#ifndef FR_SOLVER_H
#define FR_SOLVER_H

#include "euler_eq.h"
#include "type_def.h"
#include "shape_f.h"
#include "thread_pool.h"

// -----------------------------------------
// 核心类：FR 2D Euler 求解器
// -----------------------------------------

class FREulerSolver
{
public:
    FREulerSolver(Mesh&& mesh, double u_inf, double v_inf, double beta);
    FREulerSolver(Mesh&& mesh, const Vec4& init_P);
    void setGamma(double g);
    void set_fws_bc(double rho, double u, double v, double p);
    void advance(double dt);
    [[nodiscard]] double getCurrentTime() const;
    std::vector<Q9>& getNodes();
    void toVTU(const std::string& filename) const;

private:
    const Mesh mesh;
    std::vector<Q9> nodes;
    double currentTime = 0.0;
    double gamma = GAMMA;
    ThreadPool pool{std::thread::hardware_concurrency()};
    Vec4 bc = {AIR_RHO, 3 * MACH, 0.0, ONE_STD_ATM};

    void boundPreservingLimiter(std::vector<Q9>& _nodes);
    [[nodiscard]] std::vector<Q9> computeResidual(const std::vector<Q9>& _nodes);
    [[nodiscard]] Q9 computeCellResidual(const std::vector<Q9>& _nodes, int cellId) const;
    std::function<Vec4(double s)> diffFlux(const std::vector<Q9>& _nodes, int cellId, int faceType) const;
};

#endif //FR_SOLVER_H
