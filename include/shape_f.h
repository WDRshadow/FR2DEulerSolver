#ifndef SHAPE_F_H
#define SHAPE_F_H

#include "type_def.h"

#define BOTTOM 0
#define RIGHT 1
#define TOP 2
#define LEFT 3

// -----------------------------------------
// 常用数据结构及其运算符重载
// -----------------------------------------

using Q3 = std::array<Vec4, 3>; // 3个节点的Q值
using Q9 = std::array<Vec4, 9>; // 9个节点的Q值

std::vector<Q9> operator+(const std::vector<Q9>& a, const std::vector<Q9>& b);
std::vector<Q9> operator*(const std::vector<Q9>& a, double b);
void operator+=(std::vector<Q9>& a, const std::vector<Q9>& b);

// -----------------------------------------
// 3x3节点映射（逆时针）
// -----------------------------------------

std::array<int, 2> mapping_q9(int i);
std::array<int, 3> face_mapping(int faceType);
Q3 face_mapping(const Q9& q9, int faceType);

// -----------------------------------------
// 3x3 Gauss–Lobatto–Legendre 节点
// -----------------------------------------

double gll_1d(int i);
double gll_weight_1d(int i);
Vec4 gll_integrate_1d(const Q3& q3);
Point gll_2d(int i);
double gll_weight_2d(int i);
Vec4 gll_integrate_2d(const Q9& q9);

// -----------------------------------------
// 3x3 Lagrange插值函数
// -----------------------------------------

double lagrange(int i, double s);
double dlagrange(int i, double s);
Vec4 interpolate(const Q3& Qs, double s);
double lagrange(int i, double xi, double eta);
std::array<double, 2> dlagrange(int i, double xi, double eta);
Vec4 interpolate(const Q9& Qs, double xi, double eta);

// -----------------------------------------
// 2x2 双线性形函数
// -----------------------------------------

double shape(int i, double xi, double eta);
std::array<double, 2> dshape(int i, double xi, double eta);
Point interpolate(const std::array<Point, 4>& pts, double xi, double eta);

#endif //SHAPE_F_H
