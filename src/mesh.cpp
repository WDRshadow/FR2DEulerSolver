#include <cmath>

#include "mesh.h"
#include "shape_f.h"

void init_fws_mesh(Mesh& mesh, const int nx, const int ny, const double width, const double height, const double stepX,
                   const double stepH)
{
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

Matrix2d jacobian(const Mesh& mesh, int cellId, double xi, double eta)
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
