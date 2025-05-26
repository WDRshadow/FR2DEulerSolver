#include "type_def.h"

Vec4 operator+(const Vec4& a, const Vec4& b)
{
    Vec4 result;
    for (int i = 0; i < 4; ++i)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

Vec4 operator-(const Vec4& a, const Vec4& b)
{
    Vec4 result;
    for (int i = 0; i < 4; ++i)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

Vec4 operator-(const Vec4& b)
{
    Vec4 result;
    for (int i = 0; i < 4; ++i)
    {
        result[i] = -b[i];
    }
    return result;
}

Vec4 operator*(const Vec4& a, const double b)
{
    Vec4 result;
    for (int i = 0; i < 4; ++i)
    {
        result[i] = a[i] * b;
    }
    return result;
}

Vec4 operator*(const double a, const Vec4& b)
{
    return b * a;
}

Vec4 operator/(const Vec4& a, const double b)
{
    Vec4 result;
    for (int i = 0; i < 4; ++i)
    {
        result[i] = a[i] / b;
    }
    return result;
}

void operator+=(Vec4& a, const Vec4& b)
{
    for (int i = 0; i < 4; ++i)
    {
        a[i] += b[i];
    }
}

void operator-=(Vec4& a, const Vec4& b)
{
    for (int i = 0; i < 4; ++i)
    {
        a[i] -= b[i];
    }
}

Flux operator*(const Flux& a, double b)
{
    Flux result;
    for (int i = 0; i < 2; ++i)
    {
        result[i] = a[i] * b;
    }
    return result;
}

double det(const Matrix2d& J)
{
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

Matrix2d inv(const Matrix2d& J)
{
    Matrix2d J_inv;
    double detJ = det(J);
    J_inv[0][0] = J[1][1] / detJ;
    J_inv[0][1] = -J[0][1] / detJ;
    J_inv[1][0] = -J[1][0] / detJ;
    J_inv[1][1] = J[0][0] / detJ;
    return J_inv;
}

Matrix2d T(const Matrix2d& J)
{
    Matrix2d J_T;
    J_T[0][1] = J[1][0];
    J_T[1][0] = J[0][1];
    return J_T;
}

Flux operator*(const Matrix2d& A, const Flux& b)
{
    return {A[0][0] * b[0] + A[0][1] * b[1], A[1][0] * b[0] + A[1][1] * b[1]};
}
