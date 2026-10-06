// The multiprecision leg. It checks the premise that later phases build on (DESIGN D16): cpp_bin_float_50 is described
// by std::numeric_limits as an inexact, non-integer type, which is what Numerixx's default scalar_traits accept without
// an adapter. This test therefore uses Boost.Multiprecision directly, not <numerixx/adapters/multiprecision.hpp>.
// Phases 1-8 add their multiprecision instantiations here; the spike adds the root solvers (DESIGN §3.5, §7.2).
#include <numerixx/core.hpp>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <doctest/doctest.h>

#include <expected>
#include <limits>
#include <optional>
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

    TEST_CASE("cpp_bin_float_50 non-finite inputs are non_finite_input, equal ends invalid_input")
    {
        // DESIGN §6.3: the input codes do not depend on the scalar type.
        const mp50 nan = std::numeric_limits<mp50>::quiet_NaN();
        const mp50 inf = std::numeric_limits<mp50>::infinity();
        CHECK(nxx::bracket<mp50>::make(nan, mp50 { 1 }) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<mp50>::make(mp50 { 1 }, -inf) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<mp50>::make(mp50 { 1 }, mp50 { 1 }) == std::unexpected(nxx::errc::invalid_input));

        const auto res = nr::bisection {}(mp_quad, { nan, mp50 { 2 } });
        CHECK_FALSE(res.has_value());
        if (!res) {
            CHECK(res.error().code == nxx::errc::non_finite_input);
            CHECK(res.error().used == nxx::counters {});
            CHECK_FALSE(res.error().best.has_value());
        }

        const auto dx = nxx::deriv::diff(mp_quad, inf);
        CHECK_FALSE(dx.has_value());
        if (!dx) {
            CHECK(dx.error().code == nxx::errc::non_finite_input);
            CHECK(dx.error().evaluations == 0u);
        }
    }

    TEST_CASE("cpp_bin_float_50 failure estimates: the R3 order, and a pole failure carries no enclosure")
    {
        // DESIGN §6.7, §7.2: an enclosure first; the smaller width, the smaller hi/2 - lo/2 when both widths overflow;
        // then the smaller |f(x)|, with a NaN |f(x)| last.
        using est      = nr::root_estimate<mp50>;
        const mp50 nan = std::numeric_limits<mp50>::quiet_NaN();
        const mp50 inf = std::numeric_limits<mp50>::infinity();
        const mp50 m   = (std::numeric_limits<mp50>::max)();
        const auto enc = [](const mp50& lo, const mp50& hi, const mp50& fx) {
            return est { lo, fx, mp50 { hi - lo }, nr::sign_bracket<mp50> { nxx::detail::trust_me {}, lo, mp50 { -1 }, hi, mp50 { 1 } } };
        };
        const est open_one { mp50 { 0 }, mp50 { 1 } };
        const est open_nan { mp50 { 0 }, nan };
        CHECK(nxx::better_than(open_one, open_nan));
        CHECK_FALSE(nxx::better_than(open_nan, open_one));
        CHECK_FALSE(nxx::better_than(open_nan, open_nan));

        const est whole = enc(-m, m, mp50 { 0 });
        const est most  = enc(mp50 { -m / 2 }, m, mp50 { 1 });
        const est half  = enc(mp50 { 0 }, m, mp50 { 1 });
        CHECK(whole.enclosure->width() == inf);
        CHECK(most.enclosure->width() == inf);
        CHECK(nxx::better_than(most, whole));
        CHECK_FALSE(nxx::better_than(whole, most));
        CHECK(nxx::better_than(half, most));
        CHECK_FALSE(nxx::better_than(most, half));
        CHECK(nxx::better_than(enc(mp50 { 1 }, mp50 { 2 }, nan), open_one));

        const auto hyperbola = [](const mp50& x) -> mp50 { return 1 / (x - mp50 { 1 } / 3); };
        const auto res       = nr::bisection {}(hyperbola, { mp50 { 0 }, mp50 { 1 } });
        CHECK_FALSE(res.has_value());
        if (!res && res.error().best) {
            CHECK(res.error().code == nxx::errc::sign_change_not_root);
            CHECK_FALSE(res.error().best->enclosure.has_value());
            CHECK(abs(res.error().best->x - mp50 { 1 } / 3) <= res.error().best->uncertainty);
        }
        else
            FAIL_CHECK("bisection on a pole in cpp_bin_float_50 fails with a best estimate");
    }
}
