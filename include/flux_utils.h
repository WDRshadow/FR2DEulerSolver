#ifndef FLUX_UTILS_H
#define FLUX_UTILS_H

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
// Limiter
// -----------------------------------------

double minmod(double a, double b);
void bound_preserving_limiter(Q9 &cellQ, double gamma);

#endif //FLUX_UTILS_H
