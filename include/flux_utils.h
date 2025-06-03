#ifndef FLUX_UTILS_H
#define FLUX_UTILS_H

#include <Eigen/Dense>

#include "type_def.h"
#include "euler_eq.h"
#include "shape_f.h"

// -----------------------------------------
// Rusanov Numerical Flux
// -----------------------------------------

Flux rusanovFlux(const Vec4& QL, const Vec4& QR, const Point& normal, double gamma = GAMMA);

// -----------------------------------------
// Flux Gradient
// -----------------------------------------

Flux gradFlux(const Q9& q9, double xi, double eta);

// -----------------------------------------
// Flux Jacobian
// -----------------------------------------

void computeJacobianA(const std::array<double, 4>& U_avg, Eigen::Matrix4d& A, double gamma = GAMMA);
void computeJacobianB(const std::array<double, 4>& U_avg, Eigen::Matrix4d& B, double gamma = GAMMA);
Eigen::Matrix4d eigen_decompose(const Eigen::Matrix4d& A);

// -----------------------------------------
// Limiter
// -----------------------------------------

double minmod(double a, double b);
double minmod(double a, double b, double c);
void applyMinmodLimiter2D_GLL(std::vector<Q9>& _nodes, const Mesh& mesh, int cellId, double gamma);
void bound_preserving_limiter(Q9& cellQ, double gamma);

#endif //FLUX_UTILS_H
