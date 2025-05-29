find_package(GTest REQUIRED)

# -----------------------------------------

add_executable(euler_eq_test
        euler_eq_test.cpp
)
target_link_libraries(euler_eq_test
        PRIVATE
        GTest::gtest_main
        fr
)
add_test(NAME euler_eq_test COMMAND euler_eq_test)

# -----------------------------------------

add_executable(flux_utils_test
        flux_utils_test.cpp
)
target_link_libraries(flux_utils_test
        PRIVATE
        GTest::gtest_main
        fr
)
add_test(NAME flux_utils_test COMMAND flux_utils_test)

# -----------------------------------------

add_executable(mesh_test
        mesh_test.cpp
)
target_link_libraries(mesh_test
        PRIVATE
        GTest::gtest_main
        fr
)
add_test(NAME mesh_test COMMAND mesh_test)

# -----------------------------------------

add_executable(shape_f_test
        shape_f_test.cpp
)
target_link_libraries(shape_f_test
        PRIVATE
        GTest::gtest_main
        fr
)
add_test(NAME shape_f_test COMMAND shape_f_test)

# -----------------------------------------

add_executable(corr_f_test
        corr_f_test.cpp
)
target_link_libraries(corr_f_test
        PRIVATE
        GTest::gtest_main
        fr
)
add_test(NAME corr_f_test COMMAND corr_f_test)

# -----------------------------------------

add_executable(fr_solver_test
        fr_solver_test.cpp
)
target_link_libraries(fr_solver_test
    PRIVATE
    GTest::gtest_main
    fr
)
add_test(NAME fr_solver_test COMMAND fr_solver_test)
