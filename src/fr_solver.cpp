#include "fr_solver.h"
#include "corr_f.h"
#include "euler_eq.h"
#include "flux_utils.h"

using Flux3 = std::array<Flux, 3>; // 3个节点的通量

FREulerSolver::FREulerSolver(Mesh mesh, const Vec4& init_P)
    : mesh(std::move(mesh))
{
    const Vec4 Q0 = toConservative(init_P, gamma);
    nodes.resize(this->mesh.elements.size());
    for (int i = 0; i < this->mesh.elements.size(); ++i)
    {
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

void FREulerSolver::set_fws_bc(const double rho, const double u, const double p)
{
    bc = {rho, u, p};
}

void FREulerSolver::advance(const double dt)
{
    // RK4
    const auto k1 = advance(nodes, dt / 2);
    const auto k2 = advance(nodes + k1 * (dt / 2), dt / 2);
    const auto k3 = advance(nodes + k2 * (dt / 2), dt / 2);
    const auto k4 = advance(nodes + k3 * dt, dt);
    nodes += (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6);
    currentTime += dt;
}

std::vector<Q9> FREulerSolver::advance(const std::vector<Q9>& nodes, const double dt) const
{
    std::vector<Q9> newNodes = nodes;
    for (int i = 0; i < mesh.elements.size(); ++i)
    {
        auto residual = computeResidual(i);
        for (int j = 0; j < 9; ++j)
        {
            newNodes[i][j] += dt * residual[j];
        }
    }
    return newNodes;
}

void set_bc(int faceType, Q9& q9_neighbour, const Q9& q9_local, const std::array<double, 3>& bc, const double gamma)
{
    switch (faceType)
    {
    case BOTTOM:
    case TOP:
        {
            q9_neighbour = q9_local;
            q9_neighbour[0][2] = -q9_local[0][2];
            q9_neighbour[1][2] = -q9_local[1][2];
            q9_neighbour[2][2] = -q9_local[2][2];
            break;
        }
    case RIGHT:
        {
            q9_neighbour = q9_local;
            break;
        }
    case LEFT:
        {
            Vec4 bc_l = toConservative({bc[0], bc[1], 0.0, bc[2]}, gamma);
            q9_neighbour = {bc_l, bc_l, bc_l};
            break;
        }
    default:
        {
            throw std::invalid_argument("Invalid face Type");
        }
    }
}

// \hat{F}_f - F_f for one face in a cell (local coordinate system)
auto diffFlux(int cellId, const Mesh& mesh, double gamma, const std::vector<Q9>& nodes, int faceType,
              const std::array<double, 3>& bc)
{
    const auto& cell = mesh.elements[cellId];
    const auto& face = mesh.faces[cell.faceIds[faceType]];
    const auto& q9 = nodes[cellId];
    auto q9_neighbour = face.leftCell == cellId ? nodes[face.rightCell] : nodes[face.leftCell];
    if (face.leftCell == -1 || face.rightCell == -1)
    {
        set_bc(faceType, q9_neighbour, q9, bc, gamma);
    }
    auto normal = faceType == RIGHT || faceType == TOP ? face.normal : Point{-face.normal.x, -face.normal.y};
    auto points = face_mapping(faceType);
    auto points_neighbour = face_mapping((faceType + 2) % 4);
    return [=](double s)
    {
        Flux result{};
        for (int i = 0; i < 3; ++i)
        {
            auto phyFlux = physicalFlux(q9[points[i]], gamma);
            auto numFlux = rusanovFlux(q9[points[i]], q9_neighbour[points_neighbour[i]], normal, gamma);
            auto local_phyFlux = toLocalFlux(phyFlux, cellId, mesh, gll_2d(points[i]));
            auto local_numFlux = toLocalFlux(numFlux, cellId, mesh, gll_2d(points[i]));
            result[0] += (local_numFlux[0] - local_phyFlux[0]) * lagrange(i, s);
        }
        return result;
    };
}

Q9 FREulerSolver::computeResidual(const int cellId) const
{
    const auto& q9 = nodes[cellId];

    // 计算数值通量 --------------------------------------------------------
    auto dF_b = diffFlux(cellId, mesh, gamma, nodes, BOTTOM, bc);
    auto dF_r = diffFlux(cellId, mesh, gamma, nodes, RIGHT, bc);
    auto dF_t = diffFlux(cellId, mesh, gamma, nodes, TOP, bc);
    auto dF_l = diffFlux(cellId, mesh, gamma, nodes, LEFT, bc);

    // 计算重构形式 --------------------------------------------------------
    auto gradCorrFlux = [&](double xi, double eta) -> Flux
    {
        std::array<Flux, 9> fluxs{};
        for (int i = 0; i < 9; ++i)
        {
            auto F = physicalFlux(q9[i], gamma);
            Point point = gll_2d(i);
            fluxs[i] = toLocalFlux(F, cellId, mesh, point);
        }
        auto result = gradFlux(fluxs, xi, eta);
        result[0] += dF_l(eta)[0] * dg2L(xi) + dF_r(eta)[0] * dg2R(xi);
        result[1] += dF_b(xi)[1] * dg2L(eta) + dF_t(xi)[1] * dg2R(eta);
        return result;
    };

    // 计算残差 -----------------------------------------------------------
    Q9 residual;
    for (int i = 0; i < 9; ++i)
    {
        auto [xi, eta] = gll_2d(i);
        auto [dF_dxi, dG_deta] = gradCorrFlux(xi, eta);
        residual[i] = -(dF_dxi + dG_deta);
    }
    return residual;
}
