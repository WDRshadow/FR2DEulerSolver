#include <cmath>

#include "mesh.h"
#include "shape_f.h"

void init_fws_mesh(Mesh& mesh, const int nx, const int ny, const double width, const double height, const double stepX,
                   const double stepH)
{
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.width = width;
    mesh.height = height;
    auto idx = [=](int i, int j)
    {
        return j * (nx + 1) + i;
    };

    // 1. vertex
    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            double x = width * i / nx;
            double y0 = x >= stepX ? stepH : 0.0;
            double y = y0 + (height - y0) * j / ny;
            mesh.vertices.push_back({x, y});
        }
    }

    // 3. cell
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            Cell cell{};
            cell.isValid = true;
            cell.vertexIds = {idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)};
            mesh.elements.push_back(cell);
        }
    }

    // 4. face
    // 4.1 vertical plane
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            int id = j * (nx + 1) + i;
            Face face{};
            const auto& v0 = mesh.vertices[idx(i, j)];
            const auto& v1 = mesh.vertices[idx(i, j + 1)];
            double dx = v1.x - v0.x;
            double dy = v1.y - v0.y;
            double norm = std::sqrt(dx * dx + dy * dy);
            face.normal.x = dy / norm;
            face.normal.y = -dx / norm;
            if (i == 0)
            {
                face.leftCell = -1;
                face.rightCell = j * nx + 0;
            }
            else if (i == nx)
            {
                face.leftCell = j * nx + (nx - 1);
                face.rightCell = -1;
            }
            else
            {
                face.leftCell = j * nx + (i - 1);
                face.rightCell = j * nx + i;
            }
            mesh.faces.push_back(face);
            if (face.leftCell != -1)
            {
                mesh.elements[face.leftCell].faceIds[1] = id;
            }
            if (face.rightCell != -1)
            {
                mesh.elements[face.rightCell].faceIds[3] = id;
            }
        }
    }
    // 4.2 horizontal plane
    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            Face face{};
            int id = (nx + 1) * ny + nx * j + i;
            const auto& v0 = mesh.vertices[idx(i, j)];
            const auto& v1 = mesh.vertices[idx(i + 1, j)];
            double dx = v1.x - v0.x;
            double dy = v1.y - v0.y;
            double norm = std::sqrt(dx * dx + dy * dy);
            face.normal.x = -dy / norm;
            face.normal.y = dx / norm;
            if (j == 0)
            {
                face.leftCell = -1;
                face.rightCell = i;
            }
            else if (j == ny)
            {
                face.leftCell = (ny - 1) * nx + i;
                face.rightCell = -1;
            }
            else
            {
                face.leftCell = (j - 1) * nx + i;
                face.rightCell = j * nx + i;
            }
            mesh.faces.push_back(face);
            if (face.leftCell != -1)
            {
                mesh.elements[face.leftCell].faceIds[2] = id;
            }
            if (face.rightCell != -1)
            {
                mesh.elements[face.rightCell].faceIds[0] = id;
            }
        }
    }
}

void init_fws_mesh2(Mesh& mesh, const int nx, const int ny, const double width, const double height, const int stepNX,
                    const int stepNY)
{
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.width = width;
    mesh.height = height;
    mesh.vertices.clear();
    mesh.faces.clear();
    mesh.elements.clear();

    double dx = width / nx;
    double dy = height / ny;

    auto idx = [=](int i, int j)
    {
        return j * (nx + 1) + i;
    };

    // 1. 生成顶点
    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            double x = i * dx;
            double y = j * dy;
            mesh.vertices.push_back({x, y});
        }
    }

    // 2. 生成单元，跳过台阶区域
    std::vector cellMap(nx * ny, -1);
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            Cell cell{};
            int cellId = j * nx + i;
            if (i >= nx - stepNX && j < stepNY)
            {
                cell.isValid = false;
            }
            else
            {
                cell.isValid = true;
                cellMap[cellId] = cellId;
            }
            cell.vertexIds = {idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)};
            mesh.elements.push_back(cell);
        }
    }

    // 3. 生成竖直方向的边
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            int left = i == 0 ? -1 : cellMap[j * nx + (i - 1)];
            int right = i == nx ? -1 : cellMap[j * nx + i];
            if (left != -1 && right == -1 && j < stepNY)
            {
                // 如果左边有单元，右边没有单元且在台阶区域内，则标记为特殊边
                right = -2; // 特殊标记
            }
            Face face{};
            const auto& v0 = mesh.vertices[idx(i, j)];
            const auto& v1 = mesh.vertices[idx(i, j + 1)];
            double dx_ = v1.x - v0.x;
            double dy_ = v1.y - v0.y;
            double norm = std::sqrt(dx_ * dx_ + dy_ * dy_);
            face.normal.x = dy_ / norm;
            face.normal.y = -dx_ / norm;
            face.leftCell = left;
            face.rightCell = right;
            mesh.faces.push_back(face);
            int faceId = mesh.faces.size() - 1;
            if (left != -1) mesh.elements[left].faceIds[1] = faceId;
            if (right != -1) mesh.elements[right].faceIds[3] = faceId;
        }
    }

    // 4. 生成水平方向的边
    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int down = (j == 0) ? -1 : cellMap[(j - 1) * nx + i];
            int up = (j == ny) ? -1 : cellMap[j * nx + i];
            Face face{};
            const auto& v0 = mesh.vertices[idx(i, j)];
            const auto& v1 = mesh.vertices[idx(i + 1, j)];
            double dx_ = v1.x - v0.x;
            double dy_ = v1.y - v0.y;
            double norm = std::sqrt(dx_ * dx_ + dy_ * dy_);
            face.normal.x = -dy_ / norm;
            face.normal.y = dx_ / norm;
            face.leftCell = down;
            face.rightCell = up;
            mesh.faces.push_back(face);
            int faceId = mesh.faces.size() - 1;
            if (down != -1) mesh.elements[down].faceIds[2] = faceId;
            if (up != -1) mesh.elements[up].faceIds[0] = faceId;
        }
    }
}

void init_fws_mesh3(Mesh& mesh, int nx, int ny, double width, double height)
{
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.width = width;
    mesh.height = height;
    auto idx = [=](int i, int j)
    {
        return j * (nx + 1) + i;
    };

    // 1. vertex
    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            double x = width * i / nx;
            double y = height * j / ny;
            mesh.vertices.push_back({x, y});
        }
    }

    // 3. cell
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            Cell cell{};
            cell.isValid = true;
            cell.vertexIds = {idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)};
            mesh.elements.push_back(cell);
        }
    }

    // 4. face
    // 4.1 vertical plane
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            int id = j * (nx + 1) + i;
            Face face{};
            const auto& v0 = mesh.vertices[idx(i, j)];
            const auto& v1 = mesh.vertices[idx(i, j + 1)];
            double dx = v1.x - v0.x;
            double dy = v1.y - v0.y;
            double norm = std::sqrt(dx * dx + dy * dy);
            face.normal.x = dy / norm;
            face.normal.y = -dx / norm;
            if (i == 0 || i == nx)
            {
                face.leftCell = j * nx + (nx - 1);
                face.rightCell = j * nx + 0;
            }
            else
            {
                face.leftCell = j * nx + (i - 1);
                face.rightCell = j * nx + i;
            }
            mesh.faces.push_back(face);
            if (face.leftCell != -1)
            {
                mesh.elements[face.leftCell].faceIds[1] = id;
            }
            if (face.rightCell != -1)
            {
                mesh.elements[face.rightCell].faceIds[3] = id;
            }
        }
    }
    // 4.2 horizontal plane
    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            Face face{};
            int id = (nx + 1) * ny + nx * j + i;
            const auto& v0 = mesh.vertices[idx(i, j)];
            const auto& v1 = mesh.vertices[idx(i + 1, j)];
            double dx = v1.x - v0.x;
            double dy = v1.y - v0.y;
            double norm = std::sqrt(dx * dx + dy * dy);
            face.normal.x = -dy / norm;
            face.normal.y = dx / norm;
            if (j == 0 || j == ny)
            {
                face.leftCell = (ny - 1) * nx + i;
                face.rightCell = i;
            }
            else
            {
                face.leftCell = (j - 1) * nx + i;
                face.rightCell = j * nx + i;
            }
            mesh.faces.push_back(face);
            if (face.leftCell != -1)
            {
                mesh.elements[face.leftCell].faceIds[2] = id;
            }
            if (face.rightCell != -1)
            {
                mesh.elements[face.rightCell].faceIds[0] = id;
            }
        }
    }
}

Matrix2d jacobian(const Mesh& mesh, const int cellId, const double xi, const double eta)
{
    const auto& cell = mesh.elements[cellId];
    std::array<std::array<double, 2>, 2> J{};
    for (int i = 0; i < 4; ++i)
    {
        auto [dN_dxi, dN_deta] = dshape(i, xi, eta);
        J[0][0] += mesh.vertices[cell.vertexIds[i]].x * dN_dxi;
        J[0][1] += mesh.vertices[cell.vertexIds[i]].x * dN_deta;
        J[1][0] += mesh.vertices[cell.vertexIds[i]].y * dN_dxi;
        J[1][1] += mesh.vertices[cell.vertexIds[i]].y * dN_deta;
    }
    return J;
}

Matrix2d jacobian(const Mesh& mesh, const int cellId, const int faceType, const double s)
{
    const auto& cell = mesh.elements[cellId];
    double xi, eta;
    switch (faceType)
    {
    case BOTTOM:
        {
            xi = s;
            eta = -1.0;
            break;
        }
    case RIGHT:
        {
            xi = 1.0;
            eta = s;
            break;
        }
    case TOP:
        {
            xi = s;
            eta = 1.0;
            break;
        }
    case LEFT:
        {
            xi = -1.0;
            eta = s;
            break;
        }
    default:
        throw std::invalid_argument("Invalid face Type");
    }
    std::array<std::array<double, 2>, 2> J{};
    for (int i = 0; i < 4; ++i)
    {
        auto [dN_dxi, dN_deta] = dshape(i, xi, eta);
        J[0][0] += mesh.vertices[cell.vertexIds[i]].x * dN_dxi;
        J[0][1] += mesh.vertices[cell.vertexIds[i]].x * dN_deta;
        J[1][0] += mesh.vertices[cell.vertexIds[i]].y * dN_dxi;
        J[1][1] += mesh.vertices[cell.vertexIds[i]].y * dN_deta;
    }
    return J;
}
