// Phase 0: the multiprecision leg. It checks the premise that later phases build on (DESIGN D16): cpp_bin_float_50
// is described by std::numeric_limits as an inexact, non-integer type, which is what Numerixx's default
// scalar_traits (phase 1) will accept without an adapter. This test therefore uses Boost.Multiprecision directly,
// not <numerixx/adapters/multiprecision.hpp>. Phases 1-8 add their multiprecision instantiations here.
#include <numerixx/core.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <doctest/doctest.h>

#include <limits>

using mp50 = boost::multiprecision::cpp_bin_float_50;

TEST_SUITE("multiprecision")
{
    TEST_CASE("cpp_bin_float_50 has the numeric_limits the scalar traits rely on")
    {
        static_assert(std::numeric_limits<mp50>::is_specialized);
        static_assert(!std::numeric_limits<mp50>::is_integer);
        static_assert(!std::numeric_limits<mp50>::is_exact);
        static_assert(std::numeric_limits<mp50>::digits10 >= 50);
    }

    TEST_CASE("cpp_bin_float_50 arithmetic is more precise than double")
    {
        const mp50 third  = mp50 { 1 } / 3;
        const mp50 recomb = third * 3 - 1;
        CHECK(abs(recomb) < mp50 { 1e-45 });

        const mp50 root2 = sqrt(mp50 { 2 });
        CHECK(abs(root2 * root2 - 2) < mp50 { 1e-45 });
    }
}
