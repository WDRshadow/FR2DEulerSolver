#include <gtest/gtest.h>

#include "flux_utils.h"

#include "corr_f.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(NUM_FLUX, base)
{
    constexpr Vec4 P = {1.0, 1.0, 0.0, 1.0};
    const Vec4 Q = toConservative(P);
    constexpr Point normal = {1.0, 0.0};
    const Vec4 F_num = rusanovFlux(Q, Q, normal);
    const auto [F_phy, G_phy] = physicalFlux(Q);
    const auto F_n_phy = F_phy * normal.x + G_phy * normal.y;
    const auto F_diff = F_num - F_n_phy;
    EXPECT_DOUBLE_EQ(F_diff[0], 0.0);
    EXPECT_DOUBLE_EQ(F_diff[1], 0.0);
    EXPECT_DOUBLE_EQ(F_diff[2], 0.0);
    EXPECT_DOUBLE_EQ(F_diff[3], 0.0);
}

TEST(NUM_FLUX, dir)
{
    constexpr Vec4 P_local_int = {1.0, 0.0, 0.0, 1.0};
    constexpr Vec4 P_left_bc = {1.0, 1.0, 0.0, 1.0};
    const Vec4 Q_local_int = toConservative(P_local_int);
    const Vec4 Q_left_bc = toConservative(P_left_bc);
    constexpr Point normal_left = {-1.0, 0};
    auto F_phy_left = physicalFlux(Q_local_int);
    auto F_n_phy_left = F_phy_left[0] * normal_left.x + F_phy_left[1] * normal_left.y;
    auto F_num_left = rusanovFlux(Q_local_int, Q_left_bc, normal_left);
    auto F_diff_left = F_num_left - F_n_phy_left;
    EXPECT_TRUE(F_diff_left[1] < 0.0);
}

TEST(NUM_FLUX, dir2)
{
    constexpr Vec4 P_local_int = {1.0, 0.0, 0.0, 1.0};
    constexpr Vec4 P_right_bc = {1.0, 1.0, 0.0, 1.0};
    const Vec4 Q_local_int = toConservative(P_local_int);
    const Vec4 Q_right_bc = toConservative(P_right_bc);
    constexpr Point normal_right = {1.0, 0};
    auto F_phy_right = physicalFlux(Q_local_int);
    auto F_n_phy_right = F_phy_right[0] * normal_right.x + F_phy_right[1] * normal_right.y;
    auto F_num_right = rusanovFlux(Q_local_int, Q_right_bc, normal_right);
    auto F_diff_right = F_num_right - F_n_phy_right;
    EXPECT_TRUE(F_diff_right[1] < 0.0);
}

TEST(NUM_FLUX, OneD_Test)
{
    constexpr Vec4 P_local_int = {1.0, 1.0, 0.0, 1.0};
    constexpr Vec4 P_left_bc = {1.0, 1.2, 0.0, 1.0};
    const Vec4 Q_local_int = toConservative(P_local_int);
    const Vec4 Q_left_bc = toConservative(P_left_bc);
    constexpr Point normal_left = {-1.0, 0};
    // ----------------------------------------------------
    Q3 Q_local = {Q_local_int, Q_local_int, Q_local_int};
    Q3 P_local = {P_local_int, P_local_int, P_local_int};
    auto diffFlux = [&]
    {
        auto F_phy_left = physicalFlux(Q_local[0]);
        auto F_n_phy_left = F_phy_left[0] * normal_left.x + F_phy_left[1] * normal_left.y;
        auto F_num_left = rusanovFlux(Q_local[0], Q_left_bc, normal_left);
        return F_num_left - F_n_phy_left;
    };
    auto gradCorrF = [&](const double s)
    {
        Vec4 result{};
        for (int i = 0; i < 3; ++i)
        {
            auto F_phy_i = physicalFlux(Q_local[i])[0];
            result += F_phy_i * dlagrange(i, s);
        }
        result += diffFlux() * dg2R(-s);
        return result;
    };
    auto limiter = [&]
    {
        constexpr double eps = 1e-12;
        const Vec4 U_avg = 0.5 * gll_integrate_1d(Q_local);
        const auto rho_avg = U_avg[0];
        const auto p_avg = pressure(U_avg);
        double rho_min = std::numeric_limits<double>::infinity();
        double p_min = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 3; ++i)
        {
            const auto& Q = Q_local[i];
            rho_min = std::min(rho_min, Q[0]);
            p_min = std::min(p_min, pressure(Q));
        }
        double theta_rho = 1.0, theta_p = 1.0;
        if (rho_min < eps)
        {
            theta_rho = (rho_avg - eps) / (rho_avg - rho_min);
            if (theta_rho < 0.0) theta_rho = 0.0;
            if (theta_rho > 1.0) theta_rho = 1.0;
        }
        if (p_min < eps)
        {
            theta_p = (p_avg - eps) / (p_avg - p_min);
            if (theta_p < 0.0) theta_p = 0.0;
            if (theta_p > 1.0) theta_p = 1.0;
        }
        const double theta = std::min(theta_rho, theta_p);
        if (theta >= 1.0 - 1e-16)
        {
            return;
        }
        for (int i = 0; i < 3; ++i)
        {
            Vec4& U = Q_local[i];
            U = U_avg + (U - U_avg) * theta;
        }
    };
    auto update_P = [&]
    {
        for (int i = 0; i < 3; ++i)
        {
            P_local[i] = toPrimitive(Q_local[i]);
        }
    };
    auto advance = [&](const double dt)
    {
        Q3 residuals{};
        for (int i = 0; i < 3; ++i)
        {
            residuals[i] = gradCorrF(i - 1);
        }
        for (int i = 0; i < 3; ++i)
        {
            Q_local[i] -= dt * residuals[i];
        }
        limiter();
        update_P();
    };
    // ----------------------------------------------------
    double current_time = 0.0;
    for (int i = 0; i < 1000; ++i)
    {
        advance(0.1);
        current_time += 0.1;
        std::cout << "Current time: " << current_time << ". Velocity: " << P_local[0][1] << ", " << P_local[1][1] <<
            ", " << P_local[2][1] << std::endl;
        if (!(Q_local[0][0] > 0.0))
        {
            break;
        }
    }
    EXPECT_DOUBLE_EQ(P_local[0][2], 0.0);
}
