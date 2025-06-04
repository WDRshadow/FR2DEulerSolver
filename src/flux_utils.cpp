#include <cmath>

#include "flux_utils.h"
#include "shape_f.h"

double maxEig(const Vec4& Q, const Point& normal, const double gamma)
{
    auto P = toPrimitive(Q, gamma);
    double nx = normal.x;
    double ny = normal.y;
    double un = P[1] * nx + P[2] * ny;
    double c = std::sqrt(gamma * P[3] / P[0]);
    return std::abs(un) + c;
}

Vec4 rusanovFlux(const Vec4& QL, const Vec4& QR, const Point& normal, const double gamma)
{
    double sL = maxEig(QL, normal, gamma);
    double sR = maxEig(QR, normal, gamma);
    double s = std::max(sL, sR);
    auto [FL, GL] = physicalFlux(QL, gamma);
    auto [FR, GR] = physicalFlux(QR, gamma);
    double nx = normal.x, ny = normal.y;
    Vec4 Fn_nor = 0.5 * (FL * nx + GL * ny + FR * nx + GR * ny) - 0.5 * s * (QR - QL);
    return Fn_nor;
}

Flux gradFlux(const Q9& q9, const double xi, const double eta)
{
    Flux result{};
    for (int i = 0; i < 9; ++i)
    {
        auto flux = physicalFlux(q9[i]);
        auto [dl_dxi, dl_deta] = dlagrange(i, xi, eta);
        result[0] += flux[0] * dl_dxi;
        result[1] += flux[1] * dl_deta;
    }
    return result;
}

void computeJacobianA(const std::array<double, 4>& U_avg, Eigen::Matrix4d& A, const double gamma)
{
    auto [rho, rhou, rhov, E] = U_avg;
    auto [_, u, v, p] = toPrimitive(U_avg, gamma);
    const double H = (E + p) / rho;

    A.setZero();
    A(0, 0) = 0.0;
    A(0, 1) = 1.0;
    A(0, 2) = 0.0;
    A(0, 3) = 0.0;

    A(1, 0) = ((gamma - 3.0) * u * u + (gamma - 1.0) * v * v) * 0.5;
    A(1, 1) = (3.0 - gamma) * u;
    A(1, 2) = -(gamma - 1.0) * v;
    A(1, 3) = gamma - 1.0;

    A(2, 0) = (2.0 - gamma) * u * v;
    A(2, 1) = v;
    A(2, 2) = u;
    A(2, 3) = 0.0;

    A(3, 0) = u * ((gamma - 1.0) * (u * u + v * v) * 0.5 - H);
    A(3, 1) = H - (gamma - 1.0) * u * u;
    A(3, 2) = -(gamma - 1.0) * u * v;
    A(3, 3) = gamma * u;
}

void computeJacobianB(const std::array<double, 4>& U_avg, Eigen::Matrix4d& B, const double gamma)
{
    auto [rho, rhou, rhov, E] = U_avg;
    auto [_, u, v, p] = toPrimitive(U_avg, gamma);
    const double H = (E + p) / rho;

    B.setZero();
    B(0, 0) = 0.0;
    B(0, 1) = 0.0;
    B(0, 2) = 1.0;
    B(0, 3) = 0.0;

    B(1, 0) = (2.0 - gamma) * u * v;
    B(1, 1) = v;
    B(1, 2) = u;
    B(1, 3) = 0.0;

    B(2, 0) = ((gamma - 3.0) * v * v + (gamma - 1.0) * u * u) * 0.5;
    B(2, 1) = -(gamma - 1.0) * u;
    B(2, 2) = (3.0 - gamma) * v;
    B(2, 3) = gamma - 1.0;

    B(3, 0) = v * ((gamma - 1.0) * (u * u + v * v) * 0.5 - H);
    B(3, 1) = -(gamma - 1.0) * u * v;
    B(3, 2) = H - (gamma - 1.0) * v * v;
    B(3, 3) = gamma * v;
}

Eigen::Matrix4d eigen_decompose(const Eigen::Matrix4d& A)
{
    const Eigen::EigenSolver<Eigen::Matrix4d> solver(A);
    auto eigVecs = solver.eigenvectors();
    Eigen::Matrix4d R;
    for (int col = 0; col < 4; ++col)
    {
        for (int row = 0; row < 4; ++row)
        {
            R(row, col) = eigVecs(row, col).real();
        }
    }
    return R;
}

double minmod(const double a, const double b)
{
    if (a * b <= 0.0)
    {
        return 0.0;
    }
    return std::fabs(a) < std::fabs(b) ? a : b;
}

double minmod(const double a, const double b, const double c)
{
    return minmod(minmod(a, b), c);
}

void applyMinmodLimiter2D_GLL(std::vector<Q9>& _nodes, const Mesh& mesh, const int cellId, const double gamma)
{
    // Average nodes for the current cell and its neighbors
    auto& U_K_nodes = _nodes[cellId];
    const auto& U_L_nodes = _nodes[mesh.faces[mesh.elements[cellId].faceIds[LEFT]].leftCell];
    const auto& U_R_nodes = _nodes[mesh.faces[mesh.elements[cellId].faceIds[RIGHT]].rightCell];
    const auto& U_B_nodes = _nodes[mesh.faces[mesh.elements[cellId].faceIds[BOTTOM]].leftCell];
    const auto& U_T_nodes = _nodes[mesh.faces[mesh.elements[cellId].faceIds[TOP]].rightCell];
    const auto U_K_avg = gll_integrate_2d(U_K_nodes) / 4.0;
    const auto U_L_avg = gll_integrate_2d(U_L_nodes) / 4.0;
    const auto U_R_avg = gll_integrate_2d(U_R_nodes) / 4.0;
    const auto U_B_avg = gll_integrate_2d(U_B_nodes) / 4.0;
    const auto U_T_avg = gll_integrate_2d(U_T_nodes) / 4.0;

    // Slope of the cell for x and y directions
    const auto U_K_R_avg = 0.5 * gll_integrate_1d(face_mapping(U_K_nodes, RIGHT));
    const auto U_K_L_avg = 0.5 * gll_integrate_1d(face_mapping(U_K_nodes, LEFT));
    const auto U_K_T_avg = 0.5 * gll_integrate_1d(face_mapping(U_K_nodes, TOP));
    const auto U_K_B_avg = 0.5 * gll_integrate_1d(face_mapping(U_K_nodes, BOTTOM));
    const auto alpha_x = (U_K_R_avg - U_K_L_avg) / 2.0;
    const auto alpha_y = (U_K_T_avg - U_K_B_avg) / 2.0;

    // Neighborhood difference for x and y directions
    const auto delta_x_minus = U_K_avg - U_L_avg;
    const auto delta_x_plus = U_R_avg - U_K_avg;
    const auto delta_y_minus = U_K_avg - U_B_avg;
    const auto delta_y_plus = U_T_avg - U_K_avg;

    // Eigen decomposition of the Jacobian matrix A for the average state U_K_avg
    Eigen::Matrix4d A;
    computeJacobianA(U_K_avg, A, gamma);
    const auto R = eigen_decompose(A);
    const auto Rinv = R.inverse();

    // Projection to the characteristic space of A
    Vec4 beta_x_char{}, beta_y_char{};
    for (int i = 0; i < 4; ++i)
    {
        double sumbx = 0.0, sumby = 0.0;
        for (int j = 0; j < 4; ++j)
        {
            sumbx += Rinv(i, j) * alpha_x[j];
            sumby += Rinv(i, j) * alpha_y[j];
        }
        beta_x_char[i] = sumbx;
        beta_y_char[i] = sumby;
    }

    // Neighborhood differences in the characteristic space of A
    std::array<double, 4> delta_xm_char{}, delta_xp_char{}, delta_ym_char{}, delta_yp_char{};
    for (int i = 0; i < 4; ++i)
    {
        double s_xm = 0.0, s_xp = 0.0, s_ym = 0.0, s_yp = 0.0;
        for (int j = 0; j < 4; ++j)
        {
            s_xm += Rinv(i, j) * delta_x_minus[j];
            s_xp += Rinv(i, j) * delta_x_plus[j];
            s_ym += Rinv(i, j) * delta_y_minus[j];
            s_yp += Rinv(i, j) * delta_y_plus[j];
        }
        delta_xm_char[i] = s_xm;
        delta_xp_char[i] = s_xp;
        delta_ym_char[i] = s_ym;
        delta_yp_char[i] = s_yp;
    }

    // Minmod limiting in the characteristic space
    Vec4 beta_x_char_limited{}, beta_y_char_limited{};
    for (int i = 0; i < 4; ++i)
    {
        beta_x_char_limited[i] = minmod(
            beta_x_char[i],
            delta_xm_char[i],
            delta_xp_char[i]
        );
        beta_y_char_limited[i] = minmod(
            beta_y_char[i],
            delta_ym_char[i],
            delta_yp_char[i]
        );
    }

    // Projection back to the physical space of A
    Vec4 alpha_x_lim{}, alpha_y_lim{};
    for (int i = 0; i < 4; ++i)
    {
        double sum_ax = 0.0, sum_ay = 0.0;
        for (int j = 0; j < 4; ++j)
        {
            sum_ax += R(i, j) * beta_x_char_limited[j];
            sum_ay += R(i, j) * beta_y_char_limited[j];
        }
        alpha_x_lim[i] = sum_ax;
        alpha_y_lim[i] = sum_ay;
    }

    // Reconstruct the limited Q9 values
    for (int i = 0; i < 9; ++i)
    {
        const auto [xi, eta] = gll_2d(i);
        U_K_nodes[i] = U_K_avg + alpha_x_lim * xi + alpha_y_lim * eta;
    }
}


void bound_preserving_limiter(Q9& cellQ, const double gamma)
{
    constexpr double eps = 1e-12;

    // 1. rho_avg, p_avg
    const Vec4 U_avg = 0.25 * gll_integrate_2d(cellQ);
    const auto rho_avg = U_avg[0];
    const auto p_avg = pressure(U_avg, gamma);

    // 2. θ_rho, θ_p，θ = min(θ_rho, θ_p, 1.0)
    double rho_min = std::numeric_limits<double>::infinity();
    double p_min = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 9; ++i)
    {
        const auto& Q = cellQ[i];
        rho_min = std::min(rho_min, Q[0]);
        p_min = std::min(p_min, pressure(Q, gamma));
    }
    double theta_rho = 1.0, theta_p = 1.0;
    if (rho_min < eps)
    {
        theta_rho = (rho_avg - eps) / (rho_avg - rho_min);
        if (theta_rho < 0.0) theta_rho = 0.0;
        if (theta_rho > 1.0) theta_rho = 1.0;
    }
    if (p_min < eps)
    {
        theta_p = (p_avg - eps) / (p_avg - p_min);
        if (theta_p < 0.0) theta_p = 0.0;
        if (theta_p > 1.0) theta_p = 1.0;
    }
    const double theta = std::min(theta_rho, theta_p);
    if (theta >= 1.0 - 1e-16)
    {
        return;
    }

    // 3. U = U_avg + θ * (U - U_avg)
    for (int i = 0; i < 9; ++i)
    {
        Vec4& U = cellQ[i];
        U = U_avg + (U - U_avg) * theta;
    }
}
