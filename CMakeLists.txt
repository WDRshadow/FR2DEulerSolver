cmake_minimum_required(VERSION 3.31)
project(FR2DEulerSolver)

find_package(CUDAToolkit QUIET)
if(CUDAToolkit_FOUND)
    message(STATUS "CUDA detected: using GPU version.")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
else()
    message(STATUS "CUDA not detected: using CPU version.")
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

include_directories(
        PUBLIC
        include
)

add_library(fr_utils
        SHARED
        src/euler_eq.cpp
        src/type_def.cpp
        src/flux_utils.cpp
        src/mesh.cpp
        src/shape_f.cpp
        src/corr_f.cpp
)

enable_testing()
add_subdirectory(test)

add_executable(FR2DEulerSolver main.cpp)
target_link_libraries(FR2DEulerSolver fr_utils)
if(CUDAToolkit_FOUND)
    add_subdirectory(src_cuda)
    target_link_libraries(FR2DEulerSolver fr_cuda)
else()
    add_subdirectory(src_cpu)
    target_link_libraries(FR2DEulerSolver fr_cpu)
endif()
