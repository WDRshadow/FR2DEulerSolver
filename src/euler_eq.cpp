#include "euler_eq.h"

#define _rho Q[0]
#define _e Q[3]
#define _u (Q[1] / _rho)
#define _v (Q[2] / _rho)
#define _p pressure(Q, gamma)

double energy(const Vec4& P, const double gamma)
{
    return P[3] / (gamma - 1.0) + 0.5 * P[0] * (P[1] * P[1] + P[2] * P[2]);
}

double pressure(const Vec4& Q, const double gamma)
{
    return (gamma - 1.0) * (_e - 0.5 * _rho * (_u * _u + _v * _v));
}

Vec4 toConservative(const Vec4& P, const double gamma)
{
    return {P[0], P[0] * P[1], P[0] * P[2], energy(P, gamma)};
}

Vec4 toPrimitive(const Vec4& Q, const double gamma)
{
    return {_rho, _u, _v, _p};
}

Flux physicalFlux(const Vec4& Q, const double gamma)
{
    Vec4 Fx, Fy;
    Fx[0] = _rho * _u;
    Fx[1] = _rho * _u * _u + _p;
    Fx[2] = _rho * _u * _v;
    Fx[3] = _u * (_e + _p);
    Fy[0] = _rho * _v;
    Fy[1] = _rho * _u * _v;
    Fy[2] = _rho * _v * _v + _p;
    Fy[3] = _v * (_e + _p);
    return {Fx, Fy};
}
