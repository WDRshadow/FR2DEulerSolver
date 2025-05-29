cmake_minimum_required(VERSION 3.31)
project(FR2DEulerSolver)

include_directories(include)

find_package(CUDAToolkit QUIET)
if(CUDAToolkit_FOUND)
    message(STATUS "CUDA detected: using GPU version.")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    file(GLOB SRC_FILES src/*.cpp src/cuda/*.cu)
else()
    message(STATUS "CUDA not detected: using CPU version.")
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    file(GLOB SRC_FILES src/*.cpp src/cpu/*.cpp)
endif()

add_library(fr STATIC ${SRC_FILES})

add_executable(FR2DEulerSolver main.cpp)
target_link_libraries(FR2DEulerSolver fr)

enable_testing()
add_subdirectory(test)
