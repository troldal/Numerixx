// Phase 0: numerixx::linalg carries Eigen 5.0.1. The facade itself arrives in phase 5.
#include <numerixx/linalg.hpp>

// GCC 14 at -O2, and GCC 16 at -O2 without NDEBUG, report -Wnull-dereference inside Eigen 5.0.1's out-of-line
// partial_lu_impl::unblocked_lu: Block's add_to_nullable_pointer creates null branches that cannot run, because the LU
// storage is never null. -isystem does not hide middle-end warnings reported through an inline stack, but GCC checks
// the pragma state at every location in that stack, so the region must contain the first inclusion of
// PartialPivLU.h (DESIGN §5.3).
#if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wnull-dereference"
#endif
#include <Eigen/LU>
#if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic pop
#endif
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
