#ifndef FLUX_UTILS_H
#define FLUX_UTILS_H

#include "type_def.h"
#include "euler_eq.h"

// -----------------------------------------
// Rusanov Numerical Flux
// -----------------------------------------

Flux rusanovFlux(const Vec4& QL, const Vec4& QR, const Point& normal, double gamma = GAMMA);

// -----------------------------------------
// Utils
// -----------------------------------------

Flux toLocalFlux(const Flux& F_global, int cellId, const Mesh& mesh, const Point& local_point);
Flux interpolateFlux(const std::array<Flux, 9>& fluxs, double xi, double eta);
Flux gradFlux(const std::array<Flux, 9>& fluxs, double xi, double eta);

#endif //FLUX_UTILS_H
