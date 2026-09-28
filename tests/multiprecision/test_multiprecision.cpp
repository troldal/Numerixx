// The multiprecision leg. It checks the premise that later phases build on (DESIGN D16): cpp_bin_float_50 is described
// by std::numeric_limits as an inexact, non-integer type, which is what Numerixx's default scalar_traits accept without
// an adapter. This test therefore uses Boost.Multiprecision directly, not <numerixx/adapters/multiprecision.hpp>.
// Phases 1-8 add their multiprecision instantiations here; the spike adds the root solvers (DESIGN §3.5, §7.2).
#include <numerixx/core.hpp>
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <doctest/doctest.h>

#include <limits>
#include <type_traits>

using mp50 = boost::multiprecision::cpp_bin_float_50;

namespace
{
    namespace nr = nxx::roots;

    const auto mp_quad  = [](const mp50& x) -> mp50 { return x * x - 2; };
    const auto mp_dquad = [](const mp50& x) -> mp50 { return 2 * x; };

    const mp50& mp_root2()
    {
        static const mp50 value = sqrt(mp50 { 2 });
        return value;
    }

    const mp50& mp_tolerance()
    {
        static const mp50 value { "1e-45" };
        return value;
    }
}    // namespace

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

    TEST_CASE("cpp_bin_float_50 is a real scalar for the root solvers")
    {
        static_assert(nxx::real<mp50>);
        static_assert(nxx::is_real_v<mp50>);
    }

    TEST_CASE("cpp_bin_float_50 bisection converges to the square root of 2")
    {
        // floored_width{} at 168 bits needs 165 halvings of [1, 2]: the default budget (200) is sized for it, so the
        // defaults are achievable in this type too (DESIGN §3.5).
        const auto res = nr::bisection {}(mp_quad, { mp50 { 1 }, mp50 { 2 } });
        if (res) {
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(res->x)>, mp50>);
            CHECK(res->by == nr::algos::bisection);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
            CHECK(res->used.evaluations == res->used.iterations + 2u);
        }
        else
            FAIL_CHECK("bisection failed on x^2 - 2 in cpp_bin_float_50");
    }

    TEST_CASE("cpp_bin_float_50 brent converges to the square root of 2")
    {
        const auto res = nr::brent {}(mp_quad, { mp50 { 1 }, mp50 { 2 } });
        if (res) {
            CHECK(res->by == nr::algos::brent);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
        }
        else
            FAIL_CHECK("brent failed on x^2 - 2 in cpp_bin_float_50");
    }

    TEST_CASE("cpp_bin_float_50 newton converges to the square root of 2")
    {
        const auto res = nr::newton {}.with_derivative(mp_dquad)(mp_quad, mp50 { 1 });
        if (res) {
            CHECK(res->by == nr::algos::newton);
            CHECK((res->how == nxx::stop_reason::criterion || res->how == nxx::stop_reason::exact_zero));
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
            CHECK(res->used.evaluations == 2u * res->used.iterations + 1u);
        }
        else
            FAIL_CHECK("newton failed on x^2 - 2 in cpp_bin_float_50");
    }

    TEST_CASE("cpp_bin_float_50 secant converges to the square root of 2")
    {
        const auto res = nr::secant {}(mp_quad, mp50 { 1 });
        if (res) {
            CHECK(res->by == nr::algos::secant);
            CHECK((res->how == nxx::stop_reason::criterion || res->how == nxx::stop_reason::exact_zero));
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
        }
        else
            FAIL_CHECK("secant failed on x^2 - 2 in cpp_bin_float_50");
    }
}
