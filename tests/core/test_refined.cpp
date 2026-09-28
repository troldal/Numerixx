// Refined inputs (DESIGN §3.3 tiers A and B, §6.2): literals are checked at compile time by consteval constructors,
// run-time values go through make() -> std::expected, bool never converts, and roles are not interchangeable.
#include <numerixx/core.hpp>

#include <doctest/doctest.h>

#include <cstdint>
#include <expected>
#include <limits>
#include <type_traits>

namespace
{
    // Run-time values: not constant expressions, so only make() can validate them.
    double rt(double v)
    {
        volatile double x = v;
        return x;
    }

    long long rt_ll(long long v)
    {
        volatile long long x = v;
        return x;
    }

    constexpr double k_inf = std::numeric_limits<double>::infinity();
    constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();

    // A configuration parsed once into refined types (DESIGN §6.2).
    struct config
    {
        nxx::tolerance<double> tol;
        nxx::max_iterations    budget;
    };

    std::expected<config, nxx::errc> parse(double tol, long long budget)
    {
        const auto t = nxx::tolerance<double>::make(tol);
        if (!t) return std::unexpected(t.error());
        const auto b = nxx::max_iterations::make(budget);
        if (!b) return std::unexpected(b.error());
        return config { *t, *b };
    }
}    // namespace

TEST_SUITE("core")
{
    TEST_CASE("refined literals: alias-template CTAD deduces the scalar type")
    {
        constexpr nxx::tolerance t { 1e-8 };
        static_assert(std::is_same_v<decltype(t), const nxx::tolerance<double>>);
        static_assert(t.value() == 1e-8);

        constexpr nxx::tolerance tf { 1e-4f };
        static_assert(std::is_same_v<decltype(tf), const nxx::tolerance<float>>);

        constexpr nxx::tolerance tl { 1e-10L };
        static_assert(std::is_same_v<decltype(tl), const nxx::tolerance<long double>>);

        constexpr nxx::abs_tolerance a0 { 0.0 };    // >= 0: zero is legal for the absolute part of a mixed test
        static_assert(std::is_same_v<decltype(a0), const nxx::abs_tolerance<double>>);
        static_assert(a0.value() == 0.0);

        constexpr nxx::rel_tolerance r0 { 0.5 };
        static_assert(std::is_same_v<decltype(r0), const nxx::rel_tolerance<double>>);

        constexpr nxx::bracket b { 1.0, 2.0 };
        static_assert(std::is_same_v<decltype(b), const nxx::bracket<double>>);
        static_assert(b.lo() == 1.0 && b.hi() == 2.0 && b.width() == 1.0 && b.midpoint() == 1.5 && b.half_width() == 0.5);

        constexpr nxx::evaluation_budget eb { 10u };
        static_assert(eb.value() == 10u);

        CHECK(t.value() == 1e-8);
    }

    TEST_CASE("refined make(): success and failure, also in constant expressions")
    {
        // tolerance: finite and > 0
        static_assert(nxx::tolerance<double>::make(1e-8).has_value());
        static_assert(nxx::tolerance<double>::make(1e-8)->value() == 1e-8);
        static_assert(nxx::tolerance<double>::make(0.0).error() == nxx::errc::invalid_input);
        static_assert(!nxx::tolerance<double>::make(-1e-8).has_value());
        static_assert(!nxx::tolerance<double>::make(k_inf).has_value());
        static_assert(!nxx::tolerance<double>::make(k_nan).has_value());

        // abs_tolerance: finite and >= 0
        static_assert(nxx::abs_tolerance<double>::make(0.0).has_value());
        static_assert(!nxx::abs_tolerance<double>::make(-1e-300).has_value());
        static_assert(!nxx::abs_tolerance<double>::make(k_inf).has_value());
        static_assert(!nxx::abs_tolerance<double>::make(k_nan).has_value());

        // rel_tolerance: 0 <= r < 1
        static_assert(nxx::rel_tolerance<double>::make(0.0).has_value());
        static_assert(nxx::rel_tolerance<double>::make(0.999).has_value());
        static_assert(!nxx::rel_tolerance<double>::make(1.0).has_value());
        static_assert(!nxx::rel_tolerance<double>::make(-0.1).has_value());
        static_assert(!nxx::rel_tolerance<double>::make(k_inf).has_value());
        static_assert(!nxx::rel_tolerance<double>::make(k_nan).has_value());

        // evaluation_budget: >= 1
        static_assert(!nxx::evaluation_budget::make(0u).has_value());
        static_assert(nxx::evaluation_budget::make(1u).has_value());

        // The same at run time.
        const auto ok = nxx::tolerance<double>::make(rt(1e-6));
        CHECK(ok.has_value());
        if (ok) { CHECK(ok->value() == 1e-6); }
        else {
            FAIL_CHECK("tolerance<double>::make(1e-6) failed");
        }
        CHECK(nxx::tolerance<double>::make(rt(0.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::tolerance<double>::make(rt(-1.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::tolerance<double>::make(rt(k_inf)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::tolerance<double>::make(rt(k_nan)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::rel_tolerance<double>::make(rt(1.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::abs_tolerance<double>::make(rt(0.0)).has_value());

        // float and long double
        CHECK(nxx::tolerance<float>::make(1e-3f).has_value());
        CHECK_FALSE(nxx::tolerance<float>::make(std::numeric_limits<float>::infinity()).has_value());
        CHECK(nxx::tolerance<long double>::make(1e-12L).has_value());
        CHECK_FALSE(nxx::tolerance<long double>::make(-1e-12L).has_value());

        // Refined values order by their value.
        static_assert(nxx::tolerance<double> { 1e-8 } < nxx::tolerance<double> { 1e-6 });
        static_assert(nxx::tolerance<double> { 1e-8 } == *nxx::tolerance<double>::make(1e-8));
    }

    TEST_CASE("a configuration parsed once into refined types")
    {
        const auto good = parse(rt(1e-10), rt_ll(60));
        CHECK(good.has_value());
        if (good) {
            CHECK(good->tol.value() == 1e-10);
            CHECK(good->budget.value() == 60u);
        }
        else {
            FAIL_CHECK("parse(1e-10, 60) failed");
        }
        CHECK(parse(rt(-1e-10), rt_ll(60)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(parse(rt(1e-10), rt_ll(0)) == std::unexpected(nxx::errc::invalid_input));
    }

    TEST_CASE("bool never constructs a refined value")
    {
        static_assert(!std::is_constructible_v<nxx::tolerance<double>, bool>);
        static_assert(!std::is_convertible_v<bool, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::abs_tolerance<double>, bool>);
        static_assert(!std::is_constructible_v<nxx::rel_tolerance<double>, bool>);
        static_assert(!std::is_constructible_v<nxx::evaluation_budget, bool>);
        static_assert(!std::is_constructible_v<nxx::max_iterations, bool>);
        static_assert(!std::is_convertible_v<bool, nxx::max_iterations>);
        CHECK(true);
    }

    TEST_CASE("roles are not interchangeable")
    {
        static_assert(!std::is_constructible_v<nxx::tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(!std::is_convertible_v<nxx::rel_tolerance<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::rel_tolerance<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::tolerance<double>, nxx::abs_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::abs_tolerance<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::abs_tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::rel_tolerance<double>, nxx::abs_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::tolerance<double>, nxx::tolerance<float>>);

        // A criterion's single threshold is a tolerance, not a relative tolerance.
        static_assert(std::is_constructible_v<nxx::x_tol<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::x_tol<double>, nxx::rel_tolerance<double>>);
        static_assert(std::is_constructible_v<nxx::width_tol<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::rel_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::f_tol<double>, nxx::abs_tolerance<double>>);

        // An evaluation budget is not an iteration count, and vice versa.
        static_assert(std::is_constructible_v<nxx::max_evaluations, nxx::evaluation_budget>);
        static_assert(!std::is_constructible_v<nxx::max_evaluations, nxx::max_iterations>);
        static_assert(std::is_constructible_v<nxx::min_iterations, nxx::max_iterations>);
        static_assert(!std::is_constructible_v<nxx::min_iterations, nxx::evaluation_budget>);
        CHECK(true);
    }

    TEST_CASE("max_iterations: literals and make() bounds, checked in the source type")
    {
        constexpr nxx::max_iterations one { 1 };
        constexpr nxx::max_iterations top { 0xFFFF'FFFFu };
        constexpr nxx::max_iterations top_ll { 4'294'967'295LL };
        static_assert(one.value() == 1u);
        static_assert(top.value() == 4'294'967'295u);
        static_assert(top_ll == top);

        static_assert(nxx::max_iterations::make(0).error() == nxx::errc::invalid_input);
        static_assert(nxx::max_iterations::make(1)->value() == 1u);
        static_assert(nxx::max_iterations::make(2).has_value());
        static_assert(nxx::max_iterations::make(4'294'967'295LL)->value() == 4'294'967'295u);
        static_assert(!nxx::max_iterations::make(4'294'967'296LL).has_value());    // 2^32: rejected, not wrapped to 0
        static_assert(!nxx::max_iterations::make(-1).has_value());                 // rejected, not wrapped to 2^32 - 1
        static_assert(!nxx::max_iterations::make(std::numeric_limits<long long>::min()).has_value());
        static_assert(!nxx::max_iterations::make(std::numeric_limits<long long>::max()).has_value());

        CHECK(nxx::max_iterations::make(rt_ll(0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::max_iterations::make(rt_ll(1)) == nxx::max_iterations { 1 });
        CHECK(nxx::max_iterations::make(rt_ll(2)) == nxx::max_iterations { 2 });
        CHECK(nxx::max_iterations::make(rt_ll(4'294'967'295LL)) == nxx::max_iterations { 0xFFFF'FFFFu });
        CHECK(nxx::max_iterations::make(rt_ll(4'294'967'296LL)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::max_iterations::make(rt_ll(-1)) == std::unexpected(nxx::errc::invalid_input));
    }

    TEST_CASE("bracket::make re-orders, and rejects equal and non-finite endpoints")
    {
        // In constant expressions: math::isfinite takes its comparison-only path there.
        static_assert(nxx::bracket<double>::make(2.0, 1.0) == nxx::bracket { 1.0, 2.0 });
        static_assert(nxx::bracket<double>::make(1.0, 2.0) == nxx::bracket { 1.0, 2.0 });
        static_assert(nxx::bracket<double>::make(-3.0, -4.0)->lo() == -4.0);
        static_assert(nxx::bracket<double>::make(1.0, 1.0).error() == nxx::errc::invalid_input);
        static_assert(!nxx::bracket<double>::make(0.0, -0.0).has_value());    // equal values
        static_assert(!nxx::bracket<double>::make(0.0, k_inf).has_value());
        static_assert(!nxx::bracket<double>::make(-k_inf, 0.0).has_value());
        static_assert(!nxx::bracket<double>::make(k_nan, 1.0).has_value());
        static_assert(!nxx::bracket<double>::make(1.0, k_nan).has_value());
        static_assert(nxx::bracket<float>::make(2.0f, 1.0f) == nxx::bracket { 1.0f, 2.0f });
        static_assert(!nxx::bracket<float>::make(1.0f, std::numeric_limits<float>::quiet_NaN()).has_value());
        static_assert(nxx::bracket<long double>::make(2.0L, 1.0L) == nxx::bracket { 1.0L, 2.0L });

        // At run time.
        const auto b = nxx::bracket<double>::make(rt(5.0), rt(-5.0));
        CHECK(b.has_value());
        if (b) {
            CHECK(b->lo() == -5.0);
            CHECK(b->hi() == 5.0);
        }
        else {
            FAIL_CHECK("bracket<double>::make(5, -5) failed");
        }
        CHECK(nxx::bracket<double>::make(rt(1.0), rt(1.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::bracket<double>::make(rt(0.0), rt(k_inf)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::bracket<double>::make(rt(k_nan), rt(0.0)) == std::unexpected(nxx::errc::invalid_input));

        // Extreme endpoints: the midpoint and half-width stay finite.
        const auto big = nxx::bracket<double>::make(rt(-1.7e308), rt(1.7e308));
        CHECK(big.has_value());
        if (big) {
            CHECK(big->midpoint() == 0.0);
            CHECK(big->half_width() == 1.7e308);
        }
        else {
            FAIL_CHECK("bracket<double>::make(-1.7e308, 1.7e308) failed");
        }
    }

    TEST_CASE("mixed criteria: abs >= 0, 0 <= rel < 1, abs > 0 or rel > 0")
    {
        constexpr nxx::x_tol relative_only { 0.0, 1e-8 };    // purely relative: legal
        static_assert(std::is_same_v<decltype(relative_only), const nxx::x_tol<double>>);
        static_assert(relative_only.abs() == 0.0 && relative_only.rel() == 1e-8);

        constexpr nxx::x_tol single { 1e-6 };    // one threshold: a tolerance
        static_assert(std::is_same_v<decltype(single), const nxx::x_tol<double>>);
        static_assert(single.abs() == 1e-6 && single.rel() == 0.0);

        constexpr nxx::width_tol w { 0.0, 1e-3 };
        static_assert(std::is_same_v<decltype(w), const nxx::width_tol<double>>);

        static_assert(!nxx::x_tol<double>::make(0.0, 0.0).has_value());
        static_assert(nxx::x_tol<double>::make(0.0, 0.0).error() == nxx::errc::invalid_input);
        static_assert(nxx::x_tol<double>::make(0.0, 1e-8).has_value());
        static_assert(nxx::x_tol<double>::make(1e-8, 0.0).has_value());
        static_assert(!nxx::x_tol<double>::make(-1e-8, 0.1).has_value());
        static_assert(!nxx::x_tol<double>::make(1e-8, 1.0).has_value());
        static_assert(!nxx::x_tol<double>::make(k_inf, 0.1).has_value());
        static_assert(!nxx::x_tol<double>::make(k_nan, 0.1).has_value());
        static_assert(!nxx::x_tol<double>::make(0.1, k_nan).has_value());
        static_assert(!nxx::width_tol<double>::make(0.0, 0.0).has_value());
        static_assert(nxx::width_tol<double>::make(0.0, 1e-3).has_value());

        const auto xr = nxx::x_tol<double>::make(rt(0.0), rt(1e-8));
        CHECK(xr.has_value());
        if (xr) { CHECK(xr->threshold(rt(100.0)) == doctest::Approx(1e-6)); }
        else {
            FAIL_CHECK("x_tol<double>::make(0, 1e-8) failed");
        }
        CHECK(nxx::x_tol<double>::make(rt(0.0), rt(0.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::width_tol<double>::make(rt(0.0), rt(0.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::x_tol<float>::make(0.0f, 1e-4f).has_value());
        CHECK_FALSE(nxx::width_tol<long double>::make(0.0L, 0.0L).has_value());
    }

    TEST_CASE("is_refined_v: validated inputs, never plain scalars")
    {
        static_assert(nxx::is_refined_v<nxx::tolerance<double>>);
        static_assert(nxx::is_refined_v<nxx::tolerance<float>>);
        static_assert(nxx::is_refined_v<nxx::abs_tolerance<double>>);
        static_assert(nxx::is_refined_v<nxx::rel_tolerance<long double>>);
        static_assert(nxx::is_refined_v<nxx::evaluation_budget>);
        static_assert(nxx::is_refined_v<nxx::max_iterations>);
        static_assert(nxx::is_refined_v<nxx::bracket<double>>);
        static_assert(nxx::is_refined_v<const nxx::bracket<double>&>);
        static_assert(nxx::is_refined_v<nxx::x_tol<double>>);
        static_assert(nxx::is_refined_v<nxx::width_tol<double>>);
        static_assert(nxx::is_refined_v<nxx::step_tol<3, 5>>);
        static_assert(nxx::is_refined_v<nxx::floored_width>);
        static_assert(nxx::is_refined_v<nxx::f_tol<double>>);
        static_assert(nxx::is_refined_v<nxx::max_evaluations>);
        static_assert(nxx::is_refined_v<nxx::min_iterations>);
        static_assert(nxx::is_refined_v<nxx::never>);
        static_assert(nxx::is_refined_v<decltype(nxx::x_tol { 1e-6 } || nxx::max_evaluations { 10 })>);

        static_assert(!nxx::is_refined_v<double>);
        static_assert(!nxx::is_refined_v<const double&>);
        static_assert(!nxx::is_refined_v<float>);
        static_assert(!nxx::is_refined_v<int>);
        static_assert(!nxx::is_refined_v<bool>);
        static_assert(!nxx::is_refined_v<std::uint32_t>);

        // The documented pitfall (DESIGN §6.2): a consteval constructor makes is_constructible_v true although a
        // run-time double cannot construct a tolerance, so generic code must test is_refined_v instead.
        static_assert(std::is_constructible_v<nxx::tolerance<double>, double>);
        CHECK(true);
    }
}
