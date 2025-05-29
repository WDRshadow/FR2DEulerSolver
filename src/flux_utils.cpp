#include <cmath>

#include "flux_utils.h"
#include "shape_f.h"
#include "mesh.h"

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

Flux toLocalFlux(const Flux& F_global, const int cellId, const Mesh& mesh, const Point& local_point)
{
    auto J = jacobian(mesh, cellId, local_point.x, local_point.y);
    return inv(J) * F_global * det(J);
}

Flux interpolateFlux(const std::array<Flux, 9>& fluxs, const double xi, const double eta)
{
    Flux result{};
    for (int i = 0; i < 9; ++i)
    {
        result[0] += fluxs[i][0] * lagrange(i, xi, eta);
        result[1] += fluxs[i][1] * lagrange(i, xi, eta);
    }
    return result;
}

Flux gradFlux(const std::array<Flux, 9>& fluxs, double xi, double eta)
{
    Flux result{};
    for (int i = 0; i < 9; ++i)
    {
        auto [dl_dxi, dl_deta] = dlagrange(i, xi, eta);
        result[0] += fluxs[i][0] * dl_dxi;
        result[1] += fluxs[i][1] * dl_deta;
    }
    return result;
}
