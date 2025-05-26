#ifndef TYPE_DEF_H
#define TYPE_DEF_H

#include <vector>
#include <array>

#define AIR_RHO 1.225
#define ONE_STD_ATM 101325

// -----------------------------------------
// 基本数据结构
// -----------------------------------------

using Vec4 = std::array<double, 4>;
using Flux = std::array<Vec4, 2>; // Fx, Fy
using Matrix2d = std::array<std::array<double, 2>, 2>; // 2x2矩阵

struct Point
{
    double x, y;
};

struct Face
{
    int leftCell, rightCell; // 邻接单元索引, -1 表示边界
    Point normal; // 单位法向量
};

struct Cell
{
    std::array<int, 4> vertexIds; // 四边形的四个顶点
    std::array<int, 4> faceIds; // 相邻的四条边 （下，右，上，左）
};

struct Mesh
{
    std::vector<Point> vertices;
    std::vector<Face> faces;
    std::vector<Cell> elements;
};

// -----------------------------------------
// 运算符重载
// -----------------------------------------

Vec4 operator+(const Vec4& a, const Vec4& b);
Vec4 operator-(const Vec4& a, const Vec4& b);
Vec4 operator-(const Vec4& b);
Vec4 operator*(const Vec4& a, double b);
Vec4 operator*(double a, const Vec4& b);
Vec4 operator/(const Vec4& a, double b);
void operator+=(Vec4& a, const Vec4& b);
void operator-=(Vec4& a, const Vec4& b);

Flux operator*(const Flux& a, double b);

// -----------------------------------------
// 矩阵运算
// -----------------------------------------

double det(const Matrix2d& J);
Matrix2d inv(const Matrix2d& J);
Matrix2d T(const Matrix2d& J);
Flux operator*(const Matrix2d& A, const Flux& b);

#endif //TYPE_DEF_H
