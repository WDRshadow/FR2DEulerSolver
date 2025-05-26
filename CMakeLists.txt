cmake_minimum_required(VERSION 3.31)
project(FR2DEulerSolver)

set(CMAKE_CXX_STANDARD 17)

include_directories(
        include
)

add_library(fr
        STATIC
        src/euler_eq.cpp
        src/type_def.cpp
        src/flux_utils.cpp
        src/mesh.cpp
        src/shape_f.cpp
        src/corr_f.cpp
        src/fr_solver.cpp
)

add_executable(FR2DEulerSolver main.cpp)

target_link_libraries(FR2DEulerSolver fr)

enable_testing()
add_subdirectory(test)