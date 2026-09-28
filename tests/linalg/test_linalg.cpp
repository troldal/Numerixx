// Phase 0: numerixx::linalg carries Eigen 5.0.1. The facade itself arrives in phase 5.
#include <numerixx/linalg.hpp>

#include <Eigen/LU>
#include <doctest/doctest.h>

#include <cmath>

TEST_SUITE("linalg")
{
    TEST_CASE("numerixx::linalg provides Eigen 5")
    {
        // Eigen 5 keeps EIGEN_WORLD_VERSION at 3 for compatibility; its release number is EIGEN_MAJOR_VERSION.
        static_assert(EIGEN_MAJOR_VERSION == 5);

        const Eigen::Matrix2d a { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Eigen::Vector2d b { 1.0, 1.0 };
        const Eigen::Vector2d x = a.partialPivLu().solve(b);

        CHECK(std::abs(x(0) + 1.0) < 1e-12);
        CHECK(std::abs(x(1) - 1.0) < 1e-12);
    }
}
