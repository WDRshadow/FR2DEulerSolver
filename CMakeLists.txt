cmake_minimum_required(VERSION 3.31)
project(FR2DEulerSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SRC_FILES src/*.cpp)

include_directories(include)

find_package(Eigen3 REQUIRED)

add_library(fr STATIC ${SRC_FILES})
target_link_libraries(fr Eigen3::Eigen)

add_executable(FR2DEulerSolver main.cpp)
target_link_libraries(FR2DEulerSolver fr)

enable_testing()
add_subdirectory(test)
