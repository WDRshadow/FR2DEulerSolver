#ifndef MESH_H
#define MESH_H

#include "type_def.h"

// -----------------------------------------
// 前向台阶 Mesh 初始化
// -----------------------------------------

void init_fws_mesh(Mesh& mesh, int nx, int ny, double width, double height, double stepX, double stepH);
void init_fws_mesh2(Mesh& mesh, int nx, int ny, double width, double height, int stepNX, int stepNY);
void init_fws_mesh3(Mesh& mesh, int nx, int ny, double width, double height);

// -----------------------------------------
// 单元内 Jacobian 矩阵计算
// -----------------------------------------

Matrix2d jacobian(const Mesh& mesh, int cellId, double xi, double eta);
Matrix2d jacobian(const Mesh& mesh, int cellId, int faceType, double s);

#endif //MESH_H
