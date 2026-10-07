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

    // Negative requires-tests go through concepts: a negative requires-expression outside a template is ill-formed,
    // not false (DESIGN §6.2).
    template<class W, class... A>
    concept make_accepts = requires(A... a) { W::make(a...); };
    template<class... A>
    concept width_tol_ctad = requires(A... a) { nxx::width_tol { a... }; };
    template<class... A>
    concept x_tol_ctad = requires(A... a) { nxx::x_tol { a... }; };

    // Generic code that forwards a mixed tolerance (DESIGN §6.2 FLAG): through make(), never the consteval literal.
    template<class T>
    constexpr auto forward_parts(nxx::abs_tolerance<T> a, nxx::rel_tolerance<T> r)
    { return nxx::width_tol<T>::make(a, r); }
    template<class T>
    constexpr auto forward_number(T a, nxx::rel_tolerance<T> r)
    { return nxx::x_tol<T>::make(a, r); }
    // A validated tolerance with a relative part goes through the constexpr constructor, so it forwards on cl too
    // (DESIGN §6.2, §12.24).
    template<class T>
    constexpr auto forward_tolerance(nxx::tolerance<T> a, nxx::rel_tolerance<T> r)
    { return nxx::width_tol { a, r }; }
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
        // In constant expressions: math::isfinite takes its comparison-only path there. A NaN or infinite end is
        // non_finite_input, equal ends are invalid_input (DESIGN §6.3).
        static_assert(nxx::bracket<double>::make(2.0, 1.0) == nxx::bracket { 1.0, 2.0 });
        static_assert(nxx::bracket<double>::make(1.0, 2.0) == nxx::bracket { 1.0, 2.0 });
        static_assert(nxx::bracket<double>::make(-3.0, -4.0)->lo() == -4.0);
        static_assert(nxx::bracket<double>::make(1.0, 1.0).error() == nxx::errc::invalid_input);
        static_assert(!nxx::bracket<double>::make(0.0, -0.0).has_value());    // equal values
        static_assert(!nxx::bracket<double>::make(0.0, k_inf).has_value());
        static_assert(!nxx::bracket<double>::make(-k_inf, 0.0).has_value());
        static_assert(!nxx::bracket<double>::make(k_nan, 1.0).has_value());
        static_assert(!nxx::bracket<double>::make(1.0, k_nan).has_value());
        static_assert(nxx::bracket<double>::make(k_nan, 0.0).error() == nxx::errc::non_finite_input);
        static_assert(nxx::bracket<double>::make(k_inf, k_inf).error() == nxx::errc::non_finite_input);    // equal, but not finite
        static_assert(nxx::bracket<double>::make(0.0, -k_inf).error() == nxx::errc::non_finite_input);
        static_assert(nxx::bracket<double>::make(0.0, -0.0).error() == nxx::errc::invalid_input);
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
        CHECK(nxx::bracket<double>::make(rt(0.0), rt(k_inf)) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<double>::make(rt(k_nan), rt(0.0)) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<double>::make(rt(k_inf), rt(k_inf)) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<double>::make(rt(k_nan), rt(k_nan)) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<float>::make(1.0f, std::numeric_limits<float>::infinity()) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<long double>::make(std::numeric_limits<long double>::quiet_NaN(), 1.0L) ==
              std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<long double>::make(2.0L, 2.0L) == std::unexpected(nxx::errc::invalid_input));

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
        // Role-typed literals (DESIGN §6.2): one number is absolute, the relative part is always named.
        constexpr nxx::x_tol relative_only { 0.0, nxx::rel_tolerance { 1e-8 } };    // purely relative: legal
        static_assert(std::is_same_v<decltype(relative_only), const nxx::x_tol<double>>);
        static_assert(relative_only.abs() == 0.0 && relative_only.rel() == 1e-8);

        constexpr nxx::x_tol single { 1e-6 };    // one threshold: absolute
        static_assert(std::is_same_v<decltype(single), const nxx::x_tol<double>>);
        static_assert(single.abs() == 1e-6 && single.rel() == 0.0);

        constexpr nxx::width_tol w { 0.0, nxx::rel_tolerance { 1e-3 } };
        static_assert(std::is_same_v<decltype(w), const nxx::width_tol<double>>);
        constexpr nxx::width_tol mixed { 1e-10, nxx::rel_tolerance { 1e-8 } };
        static_assert(mixed.abs() == 1e-10 && mixed.rel() == 1e-8);

        // Every CTAD form: a bare absolute part, a validated one, float and long double, an integer absolute part.
        constexpr nxx::width_tol from_parts { nxx::abs_tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 } };
        static_assert(std::is_same_v<decltype(from_parts), const nxx::width_tol<double>>);
        static_assert(from_parts.abs() == mixed.abs() && from_parts.rel() == mixed.rel());
        constexpr nxx::x_tol xf { 1e-4f, nxx::rel_tolerance { 1e-3f } };
        static_assert(std::is_same_v<decltype(xf), const nxx::x_tol<float>>);
        constexpr nxx::width_tol wl { 0.0L, nxx::rel_tolerance { 1e-12L } };
        static_assert(std::is_same_v<decltype(wl), const nxx::width_tol<long double>>);
        constexpr nxx::width_tol wi { 0, nxx::rel_tolerance { 1e-8 } };
        static_assert(std::is_same_v<decltype(wi), const nxx::width_tol<double>>);

        // make(abs, rel_tolerance) and make(abs_tolerance, rel_tolerance); make(T, T) is deleted.
        static_assert(!nxx::x_tol<double>::make(0.0, nxx::rel_tolerance { 0.0 }).has_value());
        static_assert(nxx::x_tol<double>::make(0.0, nxx::rel_tolerance { 0.0 }).error() == nxx::errc::invalid_input);
        static_assert(nxx::x_tol<double>::make(0.0, nxx::rel_tolerance { 1e-8 }).has_value());
        static_assert(nxx::x_tol<double>::make(1e-8, nxx::rel_tolerance { 0.0 }).has_value());
        static_assert(!nxx::x_tol<double>::make(-1e-8, nxx::rel_tolerance { 0.1 }).has_value());
        static_assert(!nxx::x_tol<double>::make(k_inf, nxx::rel_tolerance { 0.1 }).has_value());
        static_assert(!nxx::x_tol<double>::make(k_nan, nxx::rel_tolerance { 0.1 }).has_value());
        static_assert(!nxx::width_tol<double>::make(0.0, nxx::rel_tolerance { 0.0 }).has_value());
        static_assert(nxx::width_tol<double>::make(0.0, nxx::rel_tolerance { 1e-3 }).has_value());
        static_assert(nxx::width_tol<double>::make(nxx::abs_tolerance { 0.0 }, nxx::rel_tolerance { 1e-3 }).has_value());
        static_assert(!nxx::width_tol<double>::make(nxx::abs_tolerance { 0.0 }, nxx::rel_tolerance { 0.0 }).has_value());

        // make(abs): finite and > 0, as tolerance<T>::make.
        static_assert(nxx::width_tol<double>::make(1e-10)->abs() == 1e-10);
        static_assert(nxx::width_tol<double>::make(1e-10)->rel() == 0.0);
        static_assert(nxx::x_tol<double>::make(0.0).error() == nxx::errc::invalid_input);
        static_assert(!nxx::x_tol<double>::make(-1e-10).has_value());
        static_assert(!nxx::x_tol<double>::make(k_inf).has_value());
        static_assert(!nxx::width_tol<double>::make(k_nan).has_value());

        const auto xr = nxx::rel_tolerance<double>::make(rt(1e-8)).and_then([](auto rp) { return nxx::x_tol<double>::make(rt(0.0), rp); });
        CHECK(xr.has_value());
        if (xr) { CHECK(xr->threshold(rt(100.0)) == doctest::Approx(1e-6)); }
        else {
            FAIL_CHECK("x_tol<double>::make(0, rel 1e-8) failed");
        }
        CHECK(nxx::rel_tolerance<double>::make(rt(1.0)) == std::unexpected(nxx::errc::invalid_input));    // rel >= 1
        CHECK(nxx::x_tol<double>::make(rt(0.0), nxx::rel_tolerance { 0.0 }) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::width_tol<double>::make(rt(0.0), nxx::rel_tolerance { 0.0 }) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::width_tol<double>::make(rt(0.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::x_tol<float>::make(0.0f, nxx::rel_tolerance { 1e-4f }).has_value());
        CHECK(nxx::x_tol<float>::make(1e-4f).has_value());
        CHECK_FALSE(nxx::width_tol<long double>::make(0.0L, nxx::rel_tolerance { 0.0L }).has_value());
        CHECK(nxx::width_tol<long double>::make(1e-12L).has_value());
    }

    // The run-time mirror of the literal width_tol{a, nxx::rel_tolerance{r}} checks the absolute part in-band, so a
    // mixed run-time tolerance takes two checks (DESIGN §6.2, revised on 2026-10-06, §12.21).
    TEST_CASE("mixed criteria: make(a, *rel) succeeds for a valid a and fails with invalid_input for a negative one")
    {
        const auto R = nxx::rel_tolerance<double>::make(rt(1e-8));
        CHECK(R.has_value());
        if (R) {
            const auto w = nxx::width_tol<double>::make(rt(1e-10), *R);
            CHECK(w.has_value());
            if (w) {
                CHECK(w->abs() == 1e-10);
                CHECK(w->rel() == 1e-8);
            }
            const auto x = nxx::x_tol<double>::make(rt(1e-10), *R);
            CHECK(x.has_value());
            if (x) {
                CHECK(x->abs() == 1e-10);
                CHECK(x->rel() == 1e-8);
            }
            // Purely relative at run time, and through validated role types: the same pair.
            const auto p = nxx::width_tol<double>::make(0.0, *R);
            CHECK(p.has_value());
            if (p) { CHECK(p->abs() == 0.0); }
            const auto A = nxx::abs_tolerance<double>::make(rt(1e-10));
            CHECK(A.has_value());
            if (A && w) {
                const auto v = nxx::width_tol<double>::make(*A, *R);
                CHECK(v.has_value());
                if (v) { CHECK((v->abs() == w->abs() && v->rel() == w->rel())); }
            }

            // The absolute part is checked as abs_tolerance<T>::make checks it, with the same code.
            CHECK(nxx::width_tol<double>::make(rt(-1e-10), *R) == std::unexpected(nxx::errc::invalid_input));
            CHECK(nxx::x_tol<double>::make(rt(-1e-10), *R) == std::unexpected(nxx::errc::invalid_input));
            CHECK(nxx::abs_tolerance<double>::make(rt(-1e-10)) == std::unexpected(nxx::errc::invalid_input));
            CHECK(nxx::width_tol<double>::make(rt(k_inf), *R) == std::unexpected(nxx::errc::invalid_input));
            CHECK(nxx::width_tol<double>::make(rt(k_nan), *R) == std::unexpected(nxx::errc::invalid_input));
        }
        else {
            FAIL_CHECK("rel_tolerance<double>::make(1e-8) failed");
        }
    }

    // DESIGN §6.2 FLAG: the mixed literal constructor is consteval, and cl 19.51 lacks P2564, so generic code that
    // forwards a mixed tolerance calls make(abs_tolerance, rel_tolerance) or make(T, rel_tolerance). That compiles, and
    // runs at compile time, on every compiler, cl included.
    TEST_CASE("mixed criteria: generic code forwards the parts through make")
    {
        static_assert(forward_parts(nxx::abs_tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 })->rel() == 1e-8);
        static_assert(forward_number(1e-10f, nxx::rel_tolerance { 1e-3f })->abs() == 1e-10f);
        const auto A = nxx::abs_tolerance<double>::make(rt(0.0));
        const auto R = nxx::rel_tolerance<double>::make(rt(1e-8));
        CHECK((A && R));
        if (A && R) {
            CHECK(forward_parts(*A, *R).has_value());
            CHECK(forward_number(rt(-1.0), *R) == std::unexpected(nxx::errc::invalid_input));
        }
    }

    // A validated tolerance<T> takes a relative part (DESIGN §6.2, §12.24): it is finite and > 0, so the joint
    // invariant holds without a check, and the constructor is constexpr, so run-time values and cl work too.
    TEST_CASE("mixed criteria: a validated tolerance takes a relative part, as a literal, at run time and through make")
    {
        constexpr auto w = nxx::width_tol { nxx::tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 } };
        static_assert(std::is_same_v<decltype(w), const nxx::width_tol<double>>);
        static_assert(w.abs() == 1e-10 && w.rel() == 1e-8);
        constexpr auto x = nxx::x_tol { nxx::tolerance { 1e-6f }, nxx::rel_tolerance { 1e-3f } };
        static_assert(std::is_same_v<decltype(x), const nxx::x_tol<float>>);
        static_assert(x.abs() == 1e-6f && x.rel() == 1e-3f);
        constexpr auto wl = nxx::width_tol<long double> { nxx::tolerance { 1e-12L }, nxx::rel_tolerance { 1e-9L } };
        static_assert(wl.abs() == 1e-12L && wl.rel() == 1e-9L);
        constexpr auto xe = nxx::x_tol<double> { nxx::tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 } };
        static_assert(xe.abs() == 1e-10 && xe.rel() == 1e-8);
        static_assert(nxx::width_tol<double>::make(nxx::tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 })->rel() == 1e-8);
        static_assert(nxx::x_tol<float>::make(nxx::tolerance { 1e-6f }, nxx::rel_tolerance { 1e-3f })->abs() == 1e-6f);
        static_assert(forward_tolerance(nxx::tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 }).rel() == 1e-8);
        static_assert(std::is_same_v<decltype(nxx::width_tol<double>::make(nxx::tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 })),
                                     std::expected<nxx::width_tol<double>, nxx::errc>>);

        // The other spellings keep their overloads: none of them becomes ambiguous.
        static_assert(nxx::width_tol { 1e-10, nxx::rel_tolerance { 1e-8 } }.abs() == 1e-10);
        static_assert(nxx::width_tol<double> { 0, nxx::rel_tolerance { 1e-8 } }.abs() == 0.0);
        static_assert(nxx::x_tol { nxx::abs_tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 } }.abs() == 1e-10);
        static_assert(make_accepts<nxx::width_tol<double>, nxx::tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(make_accepts<nxx::x_tol<double>, nxx::tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(make_accepts<nxx::width_tol<double>, int, nxx::rel_tolerance<double>>);
        static_assert(width_tol_ctad<nxx::tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(x_tol_ctad<nxx::tolerance<float>, nxx::rel_tolerance<float>>);

        // The roles stay typed: two tolerances, a bare relative number, swapped roles and mixed scalar types are not
        // constructible.
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::tolerance<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::tolerance<double>, double>);
        static_assert(!std::is_constructible_v<nxx::x_tol<double>, nxx::rel_tolerance<double>, nxx::tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::tolerance<float>, nxx::rel_tolerance<double>>);
        static_assert(!make_accepts<nxx::width_tol<double>, nxx::tolerance<double>, double>);
        static_assert(!make_accepts<nxx::width_tol<double>, nxx::tolerance<double>, nxx::tolerance<double>>);
        static_assert(!make_accepts<nxx::x_tol<double>, nxx::tolerance<double>, nxx::rel_tolerance<float>>);
        static_assert(!width_tol_ctad<nxx::tolerance<double>, nxx::tolerance<double>>);

        // Run-time values: the parts are validated once; the criterion then needs no check.
        const auto T = nxx::tolerance<double>::make(rt(1e-10));
        const auto R = nxx::rel_tolerance<double>::make(rt(1e-8));
        CHECK((T && R));
        if (T && R) {
            const auto a = nxx::width_tol { *T, *R };
            CHECK((a.abs() == 1e-10 && a.rel() == 1e-8));
            const auto b = nxx::x_tol<double> { *T, *R };
            CHECK((b.abs() == 1e-10 && b.rel() == 1e-8));
            const auto c = nxx::width_tol<double>::make(*T, *R);
            CHECK(c.has_value());
            if (c) { CHECK((c->abs() == 1e-10 && c->rel() == 1e-8)); }
            const auto d = nxx::x_tol<double>::make(*T, *R);
            CHECK(d.has_value());
            if (d) { CHECK((d->threshold(100.0) == b.threshold(100.0))); }
            const auto e = forward_tolerance(*T, *R);
            CHECK((e.abs() == a.abs() && e.rel() == a.rel()));
            // The same pair as the two-check path make(a, *rel).
            const auto f = nxx::width_tol<double>::make(rt(1e-10), *R);
            CHECK(f.has_value());
            if (f) { CHECK((f->abs() == a.abs() && f->rel() == a.rel())); }
        }
    }

    TEST_CASE("mixed criteria: two bare numbers and a part alone are rejected")
    {
        // Two bare numbers: deleted, with a reason (compile-fail width_tol_two_numbers, x_tol_two_numbers).
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, double, double>);
        static_assert(!std::is_constructible_v<nxx::x_tol<double>, double, double>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, double, int>);
        static_assert(!std::is_constructible_v<nxx::x_tol<float>, float, float>);
        // A part alone: deleted, with a reason (compile-fail width_tol_relative_alone).
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::rel_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::abs_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::x_tol<double>, nxx::rel_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::x_tol<double>, nxx::abs_tolerance<double>>);
        // Swapped roles, and a relative part of another scalar type: no constructor.
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, nxx::rel_tolerance<double>, double>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, double, nxx::abs_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::x_tol<double>, nxx::rel_tolerance<double>, nxx::abs_tolerance<double>>);
        static_assert(!std::is_constructible_v<nxx::width_tol<double>, double, nxx::rel_tolerance<float>>);
        // The legal mixed spellings are constructible (the constructor is consteval, so only the type is tested here).
        static_assert(std::is_constructible_v<nxx::width_tol<double>, double, nxx::rel_tolerance<double>>);
        static_assert(std::is_constructible_v<nxx::width_tol<double>, nxx::abs_tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(std::is_constructible_v<nxx::x_tol<double>, double, nxx::rel_tolerance<double>>);
        static_assert(std::is_constructible_v<nxx::x_tol<float>, float, nxx::rel_tolerance<float>>);
        // The absolute form is unchanged.
        static_assert(std::is_constructible_v<nxx::width_tol<double>, nxx::tolerance<double>>);
        static_assert(std::is_constructible_v<nxx::x_tol<double>, nxx::tolerance<double>>);

        // make: the relative part is typed (compile-fail width_tol_make_two_numbers).
        static_assert(make_accepts<nxx::width_tol<double>, double>);
        static_assert(make_accepts<nxx::width_tol<double>, double, nxx::rel_tolerance<double>>);
        static_assert(make_accepts<nxx::width_tol<double>, nxx::abs_tolerance<double>, nxx::rel_tolerance<double>>);
        static_assert(make_accepts<nxx::x_tol<double>, double, nxx::rel_tolerance<double>>);
        static_assert(!make_accepts<nxx::width_tol<double>, double, double>);
        static_assert(!make_accepts<nxx::x_tol<double>, double, double>);
        static_assert(!make_accepts<nxx::width_tol<double>, nxx::abs_tolerance<double>, double>);
        static_assert(!make_accepts<nxx::width_tol<double>, nxx::rel_tolerance<double>>);
        static_assert(!make_accepts<nxx::width_tol<double>, double, nxx::rel_tolerance<float>>);

        // Through CTAD: the deletion guides lead two numbers and a part to the deleted constructors.
        static_assert(!width_tol_ctad<double, double>);
        static_assert(!width_tol_ctad<double, int>);
        static_assert(!width_tol_ctad<nxx::rel_tolerance<double>>);
        static_assert(!width_tol_ctad<nxx::abs_tolerance<double>>);
        static_assert(!width_tol_ctad<nxx::rel_tolerance<double>, double>);
        static_assert(!x_tol_ctad<double, double>);
        static_assert(!x_tol_ctad<nxx::rel_tolerance<double>>);
        static_assert(width_tol_ctad<nxx::tolerance<double>>);
        static_assert(x_tol_ctad<nxx::tolerance<float>>);
        CHECK(true);
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
