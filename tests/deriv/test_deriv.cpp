// Numerical differentiation (DESIGN §7.1, §6.12): diff and central against a closed form, derivative_of at compile
// time, the step specifications (optimal, relative, absolute; scale 1 at x == 0) checked through the points where f is
// evaluated, fallible callbacks and their cost, a step that vanishes, cost_of, and float and long double.
#include <numerixx/deriv.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <type_traits>

namespace d = nxx::deriv;

namespace
{
    double rt(double v)
    {
        volatile double x = v;
        return x;
    }

    enum class eval_error { domain, too_big };

    // Records the points where it is evaluated (the step actually taken) and returns x^3.
    struct recorder
    {
        std::array<double, 8>* points;
        std::size_t*           count;

        double operator()(double x) const
        {
            if (*count < points->size()) (*points)[*count] = x;
            ++*count;
            return x * x * x;
        }
    };

    // The step diff takes: raw, then made exactly representable as (x + raw) - x. Separate statements, so no
    // contraction in this file either.
    double resolved(double x, double raw)
    {
        const double xh = x + raw;
        return xh - x;
    }

    // Two one-shot derivatives agree: both succeed with the same value (fault has no operator==).
    template<class R>
    bool same_value(const R& a, const R& b)
    { return a.has_value() && b.has_value() && *a == *b; }
}    // namespace

TEST_SUITE("deriv")
{
    TEST_CASE("diff and central on sin at 1 against cos(1)")
    {
        const auto   sin_fn = [](double x) { return std::sin(x); };
        const double x      = rt(1.0);
        const double exact  = std::cos(1.0);

        const auto c = d::central(sin_fn, x);
        static_assert(std::is_same_v<decltype(c), const std::expected<double, nxx::fault<nxx::none>>>);
        CHECK(c.has_value());
        CHECK(std::abs(c.value_or(0.0) - exact) < 1e-9);

        const auto c2 = d::diff(sin_fn, x);    // the defaults: central_1_2, optimal step
        CHECK(same_value(c2, c));
        const auto c12 = d::diff(sin_fn, x, d::central_1_2);
        CHECK(same_value(c12, c));

        const auto c14 = d::diff(sin_fn, x, d::central_1_4);
        CHECK(c14.has_value());
        CHECK(std::abs(c14.value_or(0.0) - exact) < 1e-11);

        // Second derivatives: -sin(1).
        const auto s22 = d::diff(sin_fn, x, d::central_2_2);
        CHECK(std::abs(s22.value_or(0.0) + std::sin(1.0)) < 1e-6);
        const auto s24 = d::diff(sin_fn, x, d::central_2_4);
        CHECK(std::abs(s24.value_or(0.0) + std::sin(1.0)) < 1e-8);

        // One-sided first order: O(h) with h = eps^(1/2).
        const auto fw = d::diff(sin_fn, x, d::forward_1_1);
        CHECK(std::abs(fw.value_or(0.0) - exact) < 1e-7);
        const auto bw = d::diff(sin_fn, x, d::backward_1_1);
        CHECK(std::abs(bw.value_or(0.0) - exact) < 1e-7);
    }

    TEST_CASE("derivative_of runs at compile time: d(x^3)/dx at 2 is 12")
    {
        constexpr auto cube = [](double x) { return x * x * x; };
        constexpr auto d3   = d::derivative_of(cube);
        constexpr auto v    = d3(2.0);
        static_assert(v.has_value());
        static_assert(nxx::math::abs(*v - 12.0) < 1e-8);
        static_assert(nxx::cost_of(d3) == 2);
        static_assert(nxx::math::abs(*d::derivative_of(cube, d::central_1_4)(2.0) - 12.0) < 1e-8);

        const auto at_run_time = d3(rt(2.0));
        CHECK(std::abs(at_run_time.value_or(0.0) - 12.0) < 1e-8);
    }

    TEST_CASE("steps: the points where f is evaluated")
    {
        std::array<double, 8> pts {};
        std::size_t           n = 0;
        const recorder        rec { &pts, &n };

        SUBCASE("optimal, central_1_2: h = 2^-17 |x|")
        {
            const double x = rt(8.0);
            CHECK(d::diff(rec, x).has_value());
            CHECK(n == 2);
            CHECK(pts[0] == 8.0 - 0x1p-14);
            CHECK(pts[1] == 8.0 + 0x1p-14);
        }
        SUBCASE("optimal, central_1_4: h = 2^-10 |x|")
        {
            CHECK(d::diff(rec, rt(1.0), d::central_1_4).has_value());
            CHECK(n == 4);
            CHECK(pts[0] == 1.0 - 0x1p-9);
            CHECK(pts[1] == 1.0 - 0x1p-10);
            CHECK(pts[2] == 1.0 + 0x1p-10);
            CHECK(pts[3] == 1.0 + 0x1p-9);
        }
        SUBCASE("optimal, central_2_2: h = 2^-13 |x|, three points")
        {
            CHECK(d::diff(rec, rt(1.0), d::central_2_2).has_value());
            CHECK(n == 3);
            CHECK(pts[0] == 1.0 - 0x1p-13);
            CHECK(pts[1] == 1.0);
            CHECK(pts[2] == 1.0 + 0x1p-13);
        }
        SUBCASE("optimal at x == 0 uses scale 1")
        {
            CHECK(d::diff(rec, rt(0.0)).has_value());
            CHECK(n == 2);
            CHECK(pts[0] == -0x1p-17);
            CHECK(pts[1] == 0x1p-17);
        }
        SUBCASE("optimal at a negative x uses |x|")
        {
            CHECK(d::diff(rec, rt(-4.0)).has_value());
            CHECK(pts[0] == -4.0 - 0x1p-15);
            CHECK(pts[1] == -4.0 + 0x1p-15);
        }
        SUBCASE("relative{f}: h = f |x|")
        {
            const double x = rt(1000.0);
            const double h = resolved(x, 1e-5 * x);
            CHECK(d::diff(rec, x, d::central_1_2, d::relative { 1e-5 }).has_value());
            CHECK(n == 2);
            CHECK(pts[0] == x - h);
            CHECK(pts[1] == x + h);
        }
        SUBCASE("relative{f} at x == 0 uses scale 1")
        {
            const double h = resolved(0.0, 1e-5);
            CHECK(d::diff(rec, rt(0.0), d::central_1_2, d::relative { 1e-5 }).has_value());
            CHECK(pts[0] == -h);
            CHECK(pts[1] == h);
        }
        SUBCASE("relative from a run-time factor")
        {
            const auto factor = nxx::tolerance<double>::make(rt(1e-4));
            CHECK(factor.has_value());
            if (factor) {
                const double x = rt(3.0);
                const double h = resolved(x, 1e-4 * x);
                CHECK(d::diff(rec, x, d::central_1_2, d::relative<double> { *factor }).has_value());
                CHECK(pts[0] == x - h);
                CHECK(pts[1] == x + h);
            }
            else {
                FAIL_CHECK("tolerance<double>::make(1e-4) failed");
            }
        }
        SUBCASE("absolute{h}: h itself, whatever x is")
        {
            const double x = rt(5.0);
            const double h = resolved(x, 1e-4);
            CHECK(d::diff(rec, x, d::central_1_2, d::absolute { 1e-4 }).has_value());
            CHECK(pts[0] == x - h);
            CHECK(pts[1] == x + h);

            n = 0;
            CHECK(d::diff(rec, rt(0.0), d::central_1_2, d::absolute { 0.25 }).has_value());
            CHECK(pts[0] == -0.25);
            CHECK(pts[1] == 0.25);
        }
    }

    TEST_CASE("steps: relative{f, typical}: h = f max(|x|, typical)")
    {
        std::array<double, 8> pts {};
        std::size_t           n = 0;
        const recorder        rec { &pts, &n };

        const double small = rt(1e-3);
        const double hs    = resolved(small, 1e-5);                                           // max(1e-3, 1) = 1
        CHECK(d::diff(rec, small, d::central_1_2, d::relative { 1e-5, 1.0 }).has_value());    // the DESIGN §7.1 spelling
        CHECK(pts[0] == small - hs);
        CHECK(pts[1] == small + hs);

        n                = 0;
        const double big = rt(100.0);
        const double hb  = resolved(big, 1e-5 * big);    // max(100, 1) = 100
        CHECK(d::diff(rec, big, d::central_1_2, d::relative { 1e-5, 1.0 }).has_value());
        CHECK(pts[0] == big - hb);
        CHECK(pts[1] == big + hb);
    }

    TEST_CASE("steps change the truncation error as expected: central_1_2 on x^3 is 3x^2 + h^2")
    {
        const auto cube = [](double x) { return x * x * x; };
        // h = 2^-4 exactly: every point and product is exact, so the result is 3 + 2^-8 exactly.
        const auto r = d::diff(cube, rt(1.0), d::central_1_2, d::absolute { 0x1p-4 });
        CHECK(r == 3.0 + 0x1p-8);
        // central_1_4 is exact on cubics.
        const auto r4 = d::diff(cube, rt(1.0), d::central_1_4, d::absolute { 0x1p-4 });
        CHECK(r4 == 3.0);
    }

    TEST_CASE("a fallible f: the fault keeps its cause and counts the evaluations so far")
    {
        std::uint32_t calls = 0;
        const auto    g     = nxx::fn::counted(
            [](double x) -> std::expected<double, eval_error> {
                if (x < 0.0) return std::unexpected(eval_error::domain);
                if (x > 1.0) return std::unexpected(eval_error::too_big);
                return x * x;
            },
            calls);
        static_assert(std::is_same_v<decltype(d::diff(g, 0.5)), std::expected<double, nxx::fault<eval_error>>>);

        SUBCASE("the first point fails")
        {
            const auto r = d::diff(g, rt(0.0), d::central_1_2, d::absolute { 0.5 });    // -0.5 fails
            CHECK_FALSE(r.has_value());
            if (!r) {
                CHECK(r.error().code == nxx::errc::callback_failed);
                CHECK(r.error().cause == eval_error::domain);
                CHECK(r.error().evaluations == 1);
                CHECK(r.error().evaluations == calls);
            }
        }
        SUBCASE("the third of four points fails")
        {
            // x = 0.8, h = 0.25: 0.3, 0.55 succeed; 1.05 fails; 1.3 is never evaluated.
            const auto r = d::diff(g, rt(0.8), d::central_1_4, d::absolute { 0.25 });
            CHECK_FALSE(r.has_value());
            if (!r) {
                CHECK(r.error().code == nxx::errc::callback_failed);
                CHECK(r.error().cause == eval_error::too_big);
                CHECK(r.error().evaluations == 3);
                CHECK(r.error().evaluations == calls);
            }
        }
        SUBCASE("success counts nothing as a fault")
        {
            const auto r = d::diff(g, rt(0.5), d::central_1_2, d::absolute { 0.25 });
            CHECK(r.has_value());
            CHECK(std::abs(r.value_or(0.0) - 1.0) < 1e-12);
            CHECK(calls == 2);
        }
        SUBCASE("through derivative_of")
        {
            const auto dg = d::derivative_of(g, d::central_1_2, d::absolute { 0.5 });
            static_assert(std::is_same_v<nxx::callback_error_t<decltype(dg), double>, eval_error>);
            const auto r = dg(rt(0.75));    // 0.25 succeeds, 1.25 fails
            CHECK_FALSE(r.has_value());
            if (!r) {
                CHECK(r.error().code == nxx::errc::callback_failed);
                CHECK(r.error().cause == eval_error::too_big);
                CHECK(r.error().evaluations == 2);
                CHECK(calls == 2);
            }
        }
    }

    TEST_CASE("a non-finite value of a plain f is non_finite_value, with its cost")
    {
        std::uint32_t calls = 0;
        const auto    f     = nxx::fn::counted([](double x) { return x > 1.0 ? std::numeric_limits<double>::quiet_NaN() : x; }, calls);
        const auto    r     = d::diff(f, rt(1.0), d::central_1_2, d::absolute { 0.5 });    // 0.5 fine, 1.5 is NaN
        CHECK_FALSE(r.has_value());
        if (!r) {
            CHECK(r.error().code == nxx::errc::non_finite_value);
            CHECK(r.error().evaluations == 2);
            CHECK(calls == 2);
        }
    }

    TEST_CASE("a step that vanishes is invalid_input, before any evaluation")
    {
        std::uint32_t calls = 0;
        const auto    f     = nxx::fn::counted([](double x) { return x; }, calls);

        const auto zero = d::diff(f, rt(1e20), d::central_1_2, d::absolute { 1e-30 });    // (1e20 + 1e-30) - 1e20 == 0
        CHECK_FALSE(zero.has_value());
        if (!zero) {
            CHECK(zero.error().code == nxx::errc::invalid_input);
            CHECK(zero.error().evaluations == 0);
        }
        CHECK(calls == 0);
    }

    TEST_CASE("a non-finite x is non_finite_input, before any evaluation")
    {
        // DESIGN §6.3, §7.1: a NaN or infinite x is an input value of its own, checked first; an overflow at a finite x
        // stays invalid_input (the cases above and below).
        constexpr double inf      = std::numeric_limits<double>::infinity();
        constexpr double nan      = std::numeric_limits<double>::quiet_NaN();
        std::uint32_t    calls    = 0;
        const auto       f        = nxx::fn::counted([](double x) { return x * x; }, calls);
        const auto       rejected = [&calls](const auto& r) {
            const bool ok = !r.has_value() && r.error().code == nxx::errc::non_finite_input && r.error().evaluations == 0 && calls == 0;
            calls         = 0;
            return ok;
        };

        CHECK(rejected(d::diff(f, rt(inf))));
        CHECK(rejected(d::diff(f, rt(-inf))));
        CHECK(rejected(d::diff(f, rt(nan))));
        CHECK(rejected(d::central(f, rt(nan))));
        CHECK(rejected(d::diff(f, rt(inf), d::forward_1_1)));                          // only x and x + h
        CHECK(rejected(d::diff(f, rt(nan), d::central_1_2, d::absolute { 1e-3 })));    // h itself would be finite
        CHECK(rejected(d::derivative_of(f)(rt(nan))));
        CHECK(rejected(d::derivative_of(f, d::central_2_4)(rt(-inf))));
        CHECK(rejected(d::numeric {}.bind(f)(rt(inf))));

        // In a constant expression, and for float and long double.
        constexpr auto sq = [](double x) { return x * x; };
        static_assert(d::diff(sq, nan).error().code == nxx::errc::non_finite_input);
        static_assert(d::diff(sq, -inf).error().code == nxx::errc::non_finite_input);
        const auto rf = d::central([](float x) { return x; }, std::numeric_limits<float>::infinity());
        CHECK_FALSE(rf.has_value());
        if (!rf) { CHECK(rf.error().code == nxx::errc::non_finite_input); }
        const auto rl = d::central([](long double x) { return x; }, std::numeric_limits<long double>::quiet_NaN());
        CHECK_FALSE(rl.has_value());
        if (!rl) { CHECK(rl.error().code == nxx::errc::non_finite_input); }

        // A fallible callback's fault type is unchanged: no cause, because f was never called.
        const auto fallible = [](double x) -> std::expected<double, int> { return x; };
        const auto rx       = d::diff(fallible, rt(nan));
        CHECK_FALSE(rx.has_value());
        if (!rx) {
            CHECK(rx.error().code == nxx::errc::non_finite_input);
            CHECK_FALSE(rx.error().cause.has_value());
        }
    }

    TEST_CASE("cost_of: the stencil's non-zero points times the cost of f")
    {
        constexpr auto sq = [](double x) { return x * x; };
        static_assert(nxx::cost_of(sq) == 1);
        static_assert(nxx::cost_of(d::derivative_of(sq)) == 2);
        static_assert(nxx::cost_of(d::derivative_of(sq, d::central_1_4)) == 4);
        static_assert(nxx::cost_of(d::derivative_of(sq, d::central_2_2)) == 3);
        static_assert(nxx::cost_of(d::derivative_of(sq, d::central_2_4)) == 5);
        static_assert(nxx::cost_of(d::derivative_of(sq, d::forward_1_1)) == 2);
        static_assert(nxx::cost_of(d::derivative_of(d::derivative_of(sq))) == 4);    // nested: 2 x 2
        static_assert(d::central_2_4.nonzero_points() == 5);

        // The derivative of a counted f costs 2 per call and makes 2 calls.
        std::uint32_t calls = 0;
        const auto    dsq   = d::derivative_of(nxx::fn::counted(sq, calls));
        CHECK(nxx::cost_of(dsq) == 2);
        CHECK(std::abs(dsq(rt(3.0)).value_or(0.0) - 6.0) < 1e-8);
        CHECK(calls == 2);
        CHECK(nxx::cost_of(std::cref(dsq)) == 2);    // through a reference_wrapper

        // The policy binds the same derivative.
        const auto bound = d::numeric {}.bind(sq);
        CHECK(nxx::cost_of(bound) == 2);
        CHECK(same_value(bound(rt(3.0)), d::derivative_of(sq)(rt(3.0))));
        CHECK(nxx::cost_of(d::numeric { d::central_1_4 }.bind(sq)) == 4);
    }

    TEST_CASE("float and long double")
    {
        const auto sinf_fn = [](float x) { return std::sin(x); };
        const auto rf      = d::central(sinf_fn, 1.0f);
        static_assert(std::is_same_v<decltype(rf), const std::expected<float, nxx::fault<nxx::none>>>);
        CHECK(std::abs(rf.value_or(0.0f) - std::cos(1.0f)) < 1e-4f);
        const auto rf4 = d::diff(sinf_fn, 1.0f, d::central_1_4);
        CHECK(std::abs(rf4.value_or(0.0f) - std::cos(1.0f)) < 1e-5f);
        const auto rf0 = d::diff(sinf_fn, 0.0f, d::central_1_2, d::relative { 1e-3f });    // scale 1 at 0
        CHECK(std::abs(rf0.value_or(0.0f) - 1.0f) < 1e-5f);

        // long double is binary64 on MSVC and x87 80-bit on MinGW: tolerances, never exact values.
        const auto sinl_fn = [](long double x) { return std::sin(x); };
        const auto rl      = d::central(sinl_fn, 1.0L);
        static_assert(std::is_same_v<decltype(rl), const std::expected<long double, nxx::fault<nxx::none>>>);
        CHECK(std::abs(rl.value_or(0.0L) - std::cos(1.0L)) < 1e-9L);
        const auto rl4 = d::diff(sinl_fn, 1.0L, d::central_1_4, d::relative { 1e-3L });
        CHECK(std::abs(rl4.value_or(0.0L) - std::cos(1.0L)) < 1e-10L);
        const auto rla = d::diff(sinl_fn, 2.0L, d::central_1_2, d::absolute { 1e-5L });
        CHECK(std::abs(rla.value_or(0.0L) - std::cos(2.0L)) < 1e-9L);

        const auto cubef = d::derivative_of([](float x) { return x * x * x; });
        CHECK(std::abs(cubef(2.0f).value_or(0.0f) - 12.0f) < 1e-2f);
        const auto cubel = d::derivative_of([](long double x) { return x * x * x; }, d::central_1_4);
        CHECK(std::abs(cubel(2.0L).value_or(0.0L) - 12.0L) < 1e-9L);
        CHECK(nxx::cost_of(cubel) == 4);
    }

    TEST_CASE("a stencil point that overflows is invalid_input, before any evaluation")
    {
        // h is finite in each call below, but a point x + k h is not. f is finite everywhere, also at +-inf (tanh(+-inf) is
        // +-1), so only diff can notice; it used to call f at +-inf and return a value that is not the stencil's estimate.
        constexpr double big      = (std::numeric_limits<double>::max)();
        std::uint32_t    calls    = 0;
        const auto       f        = nxx::fn::counted([](double x) { return std::tanh(x / (std::numeric_limits<double>::max)()); }, calls);
        const auto       rejected = [&calls](const auto& r) {
            const bool ok = !r.has_value() && r.error().code == nxx::errc::invalid_input && r.error().evaluations == 0 && calls == 0;
            calls         = 0;
            return ok;
        };

        CHECK(rejected(d::diff(f, rt(-big))));                                          // central_1_2: x - h = -inf
        CHECK(rejected(d::diff(f, rt(big))));                                           // the mirror: x + h = inf, so h is inf
        CHECK(rejected(d::diff(f, rt(big * (1 - 0x1.8p-10)), d::central_1_4)));         // x + h is finite, x + 2h = inf
        CHECK(rejected(d::diff(f, rt(0.0), d::central_1_4, d::absolute { 1e308 })));    // 2h = inf
        CHECK(rejected(d::diff(f, rt(1e308), d::central_1_4, d::relative { 0.5 })));    // x + 2h = 2e308
        CHECK(rejected(d::diff(f, rt(-big), d::central_2_4)));
        CHECK(rejected(d::diff(f, rt(-big), d::backward_1_1)));
        CHECK(rejected(d::derivative_of(f)(rt(-big))));

        // Finite points still succeed: forward_1_1 at -max samples only x and x + h.
        CHECK(d::diff(f, rt(-big), d::forward_1_1).has_value());
        CHECK(calls == 2);

        // float: x - h overflows at -max as well.
        std::uint32_t  fcalls = 0;
        const auto     ff     = nxx::fn::counted([](float x) { return std::tanh(x / (std::numeric_limits<float>::max)()); }, fcalls);
        volatile float vf     = -(std::numeric_limits<float>::max)();
        const float    xf     = vf;
        const auto     rf     = d::diff(ff, xf);
        CHECK_FALSE(rf.has_value());
        if (!rf) {
            CHECK(rf.error().code == nxx::errc::invalid_input);
            CHECK(rf.error().evaluations == 0);
        }
        CHECK(fcalls == 0);
    }
}
