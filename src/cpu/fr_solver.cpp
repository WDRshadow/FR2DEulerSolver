#include "fr_solver.h"
#include "corr_f.h"
#include "euler_eq.h"
#include "flux_utils.h"
#include "mesh.h"
#include "thread_pool.h"
#include "vtu.h"

using Flux3 = std::array<Flux, 3>;

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

void FREulerSolver::set_fws_bc(const double rho, const double u, const double p)
{
    bc = {rho, u, p};
}

void FREulerSolver::advance(const double dt)
{
    // RK4
    const auto k1 = computeResidual(nodes);
    auto u1_temp = nodes + k1 * (dt / 2);
    boundPreservingLimiter(u1_temp);
    const auto k2 = computeResidual(u1_temp);
    auto u2_temp = nodes + k2 * (dt / 2);
    boundPreservingLimiter(u2_temp);
    const auto k3 = computeResidual(u2_temp);
    auto u3_temp = nodes + k3 * dt;
    boundPreservingLimiter(u3_temp);
    const auto k4 = computeResidual(u3_temp);
    nodes += (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6);
    boundPreservingLimiter(nodes);
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
}

std::function<Flux(double s)> FREulerSolver::diffFlux(const std::vector<Q9>& _nodes, const int cellId,
                                                      const int faceType) const
{
    const auto& cell = mesh.elements[cellId];
    const auto& face = mesh.faces[cell.faceIds[faceType]];
    const auto& q9_local = _nodes[cellId];
    const auto normal = faceType == RIGHT || faceType == TOP ? face.normal : Point{-face.normal.x, -face.normal.y};
    const auto points_local = face_mapping(faceType);
    const auto points_neighbour = face_mapping((faceType + 2) % 4);
    const Q3 q3_local{
        q9_local[points_local[0]],
        q9_local[points_local[1]],
        q9_local[points_local[2]]
    };
    Q3 q3_neighbour;
    if (face.leftCell == -1 || face.rightCell == -1 || face.rightCell == -2)
    {
        switch (faceType)
        {
        case BOTTOM:
        case TOP:
            {
                q3_neighbour = q3_local;
                q3_neighbour[0][2] = -q3_neighbour[0][2];
                q3_neighbour[1][2] = -q3_neighbour[1][2];
                q3_neighbour[2][2] = -q3_neighbour[2][2];
                break;
            }
        case RIGHT:
            {
                q3_neighbour = q3_local;
                if (face.rightCell == -2)
                {
                    q3_neighbour[0][1] = -q3_neighbour[0][1];
                    q3_neighbour[1][1] = -q3_neighbour[1][1];
                    q3_neighbour[2][1] = -q3_neighbour[2][1];
                }
                break;
            }
        case LEFT:
            {
                const Vec4 bc_l = toConservative({bc[0], bc[1], 0.0, bc[2]}, gamma);
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
        q3_neighbour = {
            q9_neighbour[points_neighbour[0]],
            q9_neighbour[points_neighbour[1]],
            q9_neighbour[points_neighbour[2]]
        };
    }
    // Minmod Limiter
    // Vec4 q_avg_loc = 0.5 * gll_integrate_1d(q3_local);
    // Vec4 q_avg_neigh = 0.5 * gll_integrate_1d(q3_neighbour);
    // Q3 d_loc{};
    // Q3 d_neigh{};
    // for (int i = 0; i < 3; ++i)
    // {
    //     d_loc[i] = q3_local[i] - q_avg_loc;
    //     d_neigh[i] = q3_neighbour[i] - q_avg_neigh;
    // }
    Flux3 phyFlux{};
    Flux3 numFlux{};
    for (int i = 0; i < 3; ++i)
    {
        // Vec4 delta_loc{};
        // Vec4 delta_neigh{};
        // for (int k = 0; k < 4; ++k)
        // {
        //     delta_loc[k] = minmod(d_loc[i][k], d_loc[(i + 1) % 3][k]);
        //     delta_loc[k] = minmod(delta_loc[k], d_loc[(i + 2) % 3][k]);
        //     delta_neigh[k] = minmod(d_neigh[i][k], d_neigh[(i + 1) % 3][k]);
        //     delta_neigh[k] = minmod(delta_neigh[k], d_neigh[(i + 2) % 3][k]);
        // }
        // Vec4 QL_lim = q_avg_loc + 0.5 * delta_loc;
        // Vec4 QR_lim = q_avg_neigh - 0.5 * delta_neigh;

        phyFlux[i] = physicalFlux(q3_local[i], gamma);
        numFlux[i] = rusanovFlux(q3_local[i], q3_neighbour[i], normal, gamma);
    }
    return [=](const double s)
    {
        Flux result{};
        for (int i = 0; i < 3; ++i)
        {
            double w = lagrange(i, s);
            result[0] += (numFlux[i][0] - phyFlux[i][0]) * w;
            result[1] += (numFlux[i][1] - phyFlux[i][1]) * w;
        }
        return result;
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
        auto result = gradFlux(q9, xi, eta);
        result[0] += dF_l(eta)[0] * dg2L(xi) + dF_r(eta)[0] * dg2R(xi);
        result[1] += dF_b(xi)[1] * dg2L(eta) + dF_t(xi)[1] * dg2R(eta);
        return result;
    };

    // 计算残差 -----------------------------------------------------------
    Q9 residual;
    for (int i = 0; i < 9; ++i)
    {
        auto [xi, eta] = gll_2d(i);
        auto grad_f_corr = gradCorrFlux(xi, eta);
        const auto J = jacobian(mesh, cellId, xi, eta);
        grad_f_corr = T(inv(J)) * grad_f_corr;
        residual[i] = -(grad_f_corr[0] + grad_f_corr[1]);
    }
    return residual;
}
