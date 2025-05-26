#include "shape_f.h"

std::vector<Q9> operator+(const std::vector<Q9>& a, const std::vector<Q9>& b)
{
    if (a.size() != b.size())
    {
        throw std::invalid_argument("Vectors must be of the same size");
    }
    std::vector<Q9> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = 0; j < 9; ++j)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

std::vector<Q9> operator*(const std::vector<Q9>& a, const double b)
{
    std::vector<Q9> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = 0; j < 9; ++j)
        {
            result[i][j] = a[i][j] * b;
        }
    }
    return result;
}

void operator+=(std::vector<Q9>& a, const std::vector<Q9>& b)
{
    if (a.size() != b.size())
    {
        throw std::invalid_argument("Vectors must be of the same size");
    }
    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = 0; j < 9; ++j)
        {
            a[i][j] += b[i][j];
        }
    }
}

std::array<int, 2> mapping_q9(const int i)
{
    switch (i)
    {
    case 0:
        return {0, 0};
    case 1:
        return {1, 0};
    case 2:
        return {2, 0};
    case 3:
        return {2, 1};
    case 4:
        return {2, 2};
    case 5:
        return {1, 2};
    case 6:
        return {0, 2};
    case 7:
        return {0, 1};
    case 8:
        return {1, 1};
    default:
        throw std::invalid_argument("Invalid index for 9 node quadrilateral");
    }
}

std::array<int, 3> face_mapping(int faceType)
{
    switch (faceType)
    {
    case BOTTOM:
        return {0, 1, 2};
    case RIGHT:
        return {2, 3, 4};
    case TOP:
        return {6, 5, 4};
    case LEFT:
        return {0, 7, 6};
    default:
        throw std::invalid_argument("Invalid face Type");
    }
}

double gll_1d(const int i)
{
    switch (i)
    {
    case 0:
        return -1.0;
    case 1:
        return 0.0;
    case 2:
        return 1.0;
    default:
        throw std::invalid_argument("Invalid index for GLL node");
    }
}

Point gll_2d(const int i)
{
    auto [x, y] = mapping_q9(i);
    return {gll_1d(x), gll_1d(y)};
}

double lagrange(const int i, const double s)
{
    switch (i)
    {
    case 0:
        return s * (s - 1.0) / 2.0;
    case 1:
        return 1.0 - s * s;
    case 2:
        return s * (s + 1.0) / 2.0;
    default:
        throw std::invalid_argument("Invalid index for Lagrange interpolation");
    }
}

double dlagrange(const int i, const double s)
{
    switch (i)
    {
    case 0:
        return s - 0.5;
    case 1:
        return -2.0 * s;
    case 2:
        return s + 0.5;
    default:
        throw std::invalid_argument("Invalid index for Lagrange interpolation");
    }
}

Vec4 interpolate(const Q3& Qs, const double s)
{
    Vec4 result;
    for (int i = 0; i < 3; ++i)
    {
        result += Qs[i] * lagrange(i, s);
    }
    return result;
}

double lagrange(const int i, const double xi, const double eta)
{
    auto [_i, _j] = mapping_q9(i);
    return lagrange(_i, xi) * lagrange(_j, eta);
}

std::array<double, 2> dlagrange(const int i, const double xi, const double eta)
{
    auto [_i, _j] = mapping_q9(i);
    return {
        dlagrange(_i, xi) * lagrange(_j, eta),
        lagrange(_i, xi) * dlagrange(_j, eta)
    };
}

Vec4 interpolate(const Q9& Qs, const double xi, const double eta)
{
    Vec4 result;
    for (int i = 0; i < 9; ++i)
    {
        result += Qs[i] * lagrange(i, xi, eta);
    }
    return result;
}

double shape(const int i, const double xi, const double eta)
{
    switch (i)
    {
    case 0:
        return 1.0 / 4.0 * (1 - xi) * (1 - eta);
    case 1:
        return 1.0 / 4.0 * (1 + xi) * (1 - eta);
    case 2:
        return 1.0 / 4.0 * (1 + xi) * (1 + eta);
    case 3:
        return 1.0 / 4.0 * (1 - xi) * (1 + eta);
    default:
        throw std::invalid_argument("Invalid index for shape function");
    }
}

std::array<double, 2> dshape(const int i, const double xi, const double eta)
{
    switch (i)
    {
    case 0:
        return {-1.0 / 4.0 * (1 - eta), -1.0 / 4.0 * (1 - xi)};
    case 1:
        return {1.0 / 4.0 * (1 - eta), -1.0 / 4.0 * (1 + xi)};
    case 2:
        return {1.0 / 4.0 * (1 + eta), 1.0 / 4.0 * (1 + xi)};
    case 3:
        return {-1.0 / 4.0 * (1 + eta), 1.0 / 4.0 * (1 - xi)};
    default:
        throw std::invalid_argument("Invalid index for shape function");
    }
}
