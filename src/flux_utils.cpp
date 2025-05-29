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

Flux rusanovFlux(const Vec4& QL, const Vec4& QR, const Point& normal, const double gamma)
{
    double sL = maxEig(QL, normal, gamma);
    double sR = maxEig(QR, normal, gamma);
    double s = std::max(sL, sR);
    auto [FL, GL] = physicalFlux(QL, gamma);
    auto [FR, GR] = physicalFlux(QR, gamma);
    double nx = normal.x, ny = normal.y;
    Vec4 Fn_nor = 0.5 * (FL * nx + GL * ny + FR * nx + GR * ny) - 0.5 * s * (QR - QL);
    Flux Fn = {Fn_nor * nx, Fn_nor * ny};
    return Fn;
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

double minmod(const double a, const double b)
{
    if (a * b <= 0.0)
    {
        return 0.0;
    }
    return std::fabs(a) < std::fabs(b) ? a : b;
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
