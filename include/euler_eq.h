#ifndef EULER_EQ_H
#define EULER_EQ_H

#include "type_def.h"

#define GAMMA 1.4

// -----------------------------------------
// 能量/压力计算
// -----------------------------------------

double energy(const Vec4& P, double gamma = GAMMA);
double pressure(const Vec4& Q, double gamma = GAMMA);

// -----------------------------------------
// 保守量/原始量转换
// -----------------------------------------

Vec4 toPrimitive(const Vec4& Q, double gamma = GAMMA);
Vec4 toConservative(const Vec4& P, double gamma = GAMMA);

// -----------------------------------------
// 物理通量计算
// -----------------------------------------

Flux physicalFlux(const Vec4& Q, double gamma = GAMMA);

#endif //EULER_EQ_H
