# Flux Reconstruction 2D Euler Solver

## Algorithm

1. **Mesh Generation**: quadrilateral，Gauss-Lobatto-Legendre (3x3)
2. **Numerical Flux**: Rusanov
3. **Correction Function**: g2
4. **Time Integration**: RK4

## Structure

```text
FR2DEulerSolver                                 // 2D Euler 求解器
├── gamma: double                               // 比热比
├── mesh: Mesh                                  // 网格
    ├── vertices: std::vector<Point>            // 网格几何顶点
        ├── x: double                           // 点x坐标
        └── y: double                           // 点y坐标
    ├── faces: std::vector<Face>                // 通量边界、边界类型、法向量
        ├── normal: Point                       // 法向量
        ├── leftCell: int                       // 左侧单元index
        └── rightCell: int                      // 右侧单元index
    └── elements: std::vector<Cell>             // 控制体，连接到节点和边
        ├── vertexIds: std::array<int, 4>       // 网格单元几何顶点index
        └── faceIds: std::array<int, 4>         // 网格单元几何边index
└── nodes: std::vector<Q9>                      // 各个单元的真实自由度（Q 点）
    └── Q: std::array<Vec4, 9>                  // 点x坐标
        ├── rho: double                         // 密度
        ├── rou*u: double                       // x 方向动量密度
        ├── rou*v: double                       // y 方向动量密度
        └── E: double                           // 能量
```
