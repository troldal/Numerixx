// Phase 0: numerixx::multiprecision provides the adapter header and Boost.Multiprecision. The adapter's content
// arrives in phase 9 (DESIGN §10.3).
#include <numerixx/adapters/multiprecision.hpp>

#include <doctest/doctest.h>

TEST_SUITE("multiprecision")
{
    TEST_CASE("the multiprecision adapter header provides cpp_bin_float_50")
    {
        using mp50         = boost::multiprecision::cpp_bin_float_50;
        const mp50 residue = mp50 { 1 } / 3 * 3 - 1;
        CHECK(abs(residue) < mp50 { 1e-45 });
    }
}
