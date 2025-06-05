#include "fr_solver.h"
#include "corr_f.h"
#include "euler_eq.h"
#include "flux_utils.h"
#include "mesh.h"
#include "thread_pool.h"
#include "vtu.h"

using Flux3 = std::array<Flux, 3>;

FREulerSolver::FREulerSolver(Mesh&& mesh, const double u_inf, const double v_inf, const double beta)
    : mesh(std::move(mesh))
{
    const double cx = mesh.width / 2;
    const double cy = mesh.height / 2;
    auto getP = [&](int cellId, int nodeId) -> Vec4
    {
        auto [xi, eta] = gll_2d(nodeId);
        auto [x, y] = interpolate({
                                      this->mesh.vertices[this->mesh.elements[cellId].vertexIds[0]],
                                      this->mesh.vertices[this->mesh.elements[cellId].vertexIds[1]],
                                      this->mesh.vertices[this->mesh.elements[cellId].vertexIds[2]],
                                      this->mesh.vertices[this->mesh.elements[cellId].vertexIds[3]]
                                  }, xi, eta);
        double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
        double T = 1 - (gamma - 1) * beta * beta / (8 * gamma * EIGEN_PI * EIGEN_PI) * std::exp(1 - r * r);
        double rho = std::pow(T, 1.0 / (gamma - 1.0));
        double u = u_inf - beta / (2 * EIGEN_PI) * std::exp((1 - r * r) / 2) * (y - cy);
        double v = v_inf + beta / (2 * EIGEN_PI) * std::exp((1 - r * r) / 2) * (x - cx);
        double p = std::pow(rho, gamma);

        return {rho, u, v, p};
    };
    nodes.resize(this->mesh.elements.size());
    for (int i = 0; i < this->mesh.elements.size(); ++i)
    {
        if (!this->mesh.elements[i].isValid) continue;
        for (int j = 0; j < 9; ++j)
        {
            nodes[i][j] = toConservative(getP(i, j), gamma);
        }
    }
}

FREulerSolver::FREulerSolver(Mesh&& mesh, const Vec4& init_P)
    : mesh(std::move(mesh))
{
    const Vec4 Q0 = toConservative(init_P, gamma);
    nodes.resize(this->mesh.elements.size());
    for (int i = 0; i < this->mesh.elements.size(); ++i)
    {
        if (!this->mesh.elements[i].isValid) continue;
        for (int j = 0; j < 9; ++j)
        {
            nodes[i][j] = Q0;
        }
    }
}

void FREulerSolver::setGamma(const double g)
{
    gamma = g;
}

void FREulerSolver::set_fws_bc(const double rho, const double u, const double v, const double p)
{
    bc = {rho, u, v, p};
}

void FREulerSolver::advance(const double dt)
{
    // RK4
    const auto k1 = computeResidual(nodes);
    auto u1_temp = nodes + k1 * (dt / 2);
    // boundPreservingLimiter(u1_temp);
    const auto k2 = computeResidual(u1_temp);
    auto u2_temp = nodes + k2 * (dt / 2);
    // boundPreservingLimiter(u2_temp);
    const auto k3 = computeResidual(u2_temp);
    auto u3_temp = nodes + k3 * dt;
    // boundPreservingLimiter(u3_temp);
    const auto k4 = computeResidual(u3_temp);
    nodes += (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6);
    // boundPreservingLimiter(nodes);

    // Euler
    // const auto residual = computeResidual(nodes);
    // nodes += residual * dt;
    // boundPreservingLimiter(nodes);

    currentTime += dt;
}

double FREulerSolver::getCurrentTime() const
{
    return currentTime;
}

std::vector<Q9> FREulerSolver::getNodes() const
{
    return nodes;
}

void FREulerSolver::boundPreservingLimiter(std::vector<Q9>& _nodes)
{
    std::vector<std::function<void()>> tasks;
    tasks.reserve(mesh.elements.size());
    for (int i = 0; i < mesh.elements.size(); ++i)
    {
        if (mesh.elements[i].isValid == false) continue;
        tasks.emplace_back([&, i]
        {
            bound_preserving_limiter(_nodes[i], gamma);
        });
    }
    pool.enqueue_bulk(tasks);
    pool.join();
}

std::vector<Q9> FREulerSolver::computeResidual(const std::vector<Q9>& _nodes)
{
    std::vector<Q9> k(nodes.size());
    std::vector<std::function<void()>> tasks;
    tasks.reserve(mesh.elements.size());
    for (int i = 0; i < mesh.elements.size(); ++i)
    {
        if (mesh.elements[i].isValid == false) continue;
        tasks.emplace_back([&, i]
        {
            auto residual = computeCellResidual(_nodes, i);
            k[i] = residual;
        });
    }
    pool.enqueue_bulk(tasks);
    pool.join();
    return k;
}

void FREulerSolver::toVTU(const std::string& filename) const
{
    writeFRSolutionVTU(filename, nodes, mesh, gamma);
    // writeAvgFRSolutionVTU(filename, nodes, mesh, gamma);
}

std::function<Vec4(double s)> FREulerSolver::diffFlux(const std::vector<Q9>& _nodes, const int cellId,
                                                      const int faceType) const
{
    const auto& cell = mesh.elements[cellId];
    const auto& face = mesh.faces[cell.faceIds[faceType]];
    const auto& q9_local = _nodes[cellId];
    const auto normal = faceType == RIGHT || faceType == TOP ? face.normal : Point{-face.normal.x, -face.normal.y};
    const Q3 q3_local = face_mapping(q9_local, faceType);
    Q3 q3_neighbour;
    if (face.leftCell < 0 || face.rightCell < 0)
    {
        auto boundaryType = face.leftCell < 0 ? face.leftCell : face.rightCell;
        switch (boundaryType)
        {
        case X_WALL:
            {
                q3_neighbour = q3_local;
                q3_neighbour[0][1] = -q3_neighbour[0][1];
                q3_neighbour[1][1] = -q3_neighbour[1][1];
                q3_neighbour[2][1] = -q3_neighbour[2][1];
            }
        case Y_WALL:
            {
                q3_neighbour = q3_local;
                q3_neighbour[0][2] = -q3_neighbour[0][2];
                q3_neighbour[1][2] = -q3_neighbour[1][2];
                q3_neighbour[2][2] = -q3_neighbour[2][2];
                break;
            }
        case OUTLET:
            {
                q3_neighbour = q3_local;
                break;
            }
        case INLET:
            {
                const Vec4 bc_l = toConservative(bc, gamma);
                q3_neighbour = {bc_l, bc_l, bc_l};
                break;
            }
        default:
            {
            }
        }
    }
    else
    {
        const auto& q9_neighbour = face.leftCell == cellId ? nodes[face.rightCell] : nodes[face.leftCell];
        q3_neighbour = face_mapping(q9_neighbour, (faceType + 2) % 4);
    }
    Q3 phyFlux{};
    Q3 numFlux{};
    for (int i = 0; i < 3; ++i)
    {
        auto physicalFlux_i = physicalFlux(q3_local[i], gamma);
        phyFlux[i] = physicalFlux_i[0] * normal.x + physicalFlux_i[1] * normal.y;
        numFlux[i] = rusanovFlux(q3_local[i], q3_neighbour[i], normal, gamma);
    }
    const auto Js = jacobian(mesh, cellId, faceType);
    return [=](const double s)
    {
        const int i = static_cast<int>(s + 1);
        return (numFlux[i] - phyFlux[i]) * Js;
    };
}

Q9 FREulerSolver::computeCellResidual(const std::vector<Q9>& _nodes, const int cellId) const
{
    // 计算数值通量 --------------------------------------------------------
    const auto dF_b = diffFlux(_nodes, cellId, BOTTOM);
    const auto dF_r = diffFlux(_nodes, cellId, RIGHT);
    const auto dF_t = diffFlux(_nodes, cellId, TOP);
    const auto dF_l = diffFlux(_nodes, cellId, LEFT);

    // 计算重构形式 --------------------------------------------------------
    const auto& q9 = _nodes[cellId];
    auto gradCorrFlux = [&](const double xi, const double eta)
    {
        const auto J = jacobian(mesh, cellId, xi, eta);
        const auto J_det = det(J);
        auto result = gradFlux(q9, xi, eta);
        result[0] += (dF_l(eta) * dg2R(-xi) + dF_r(eta) * dg2R(xi)) / J_det;
        result[1] += (dF_b(xi) * dg2R(-eta) + dF_t(xi) * dg2R(eta)) / J_det;
        return T(inv(J)) * result;
    };

    // 计算残差 -----------------------------------------------------------
    Q9 residual;
    for (int i = 0; i < 9; ++i)
    {
        auto [xi, eta] = gll_2d(i);
        auto grad_f_corr = gradCorrFlux(xi, eta);
        residual[i] = -(grad_f_corr[0] + grad_f_corr[1]);
    }
    return residual;
}
