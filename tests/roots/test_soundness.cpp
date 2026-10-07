// Criterion soundness (DESIGN §9.3; spike exit criterion 7, §10.2): on every spike solver run standalone, a success
// with stop_reason::criterion implies the criterion's guarantee for the returned estimate. Also: a bracketing solver
// never reports convergence while its width exceeds the requested tolerance; resolution_limit implies adjacent
// endpoints; every intermediate bracket keeps a sign change and strictly shrinks; x_tol on a bracketing solver does
// not compile.
//
// The random problems come from a fixed-seed std::mt19937. The standard libraries draw different sequences from the
// distributions, so the problems, and the number of stops of each kind, differ between toolchains; the guarantees
// must hold on all of them. The minimum counts only make sure that a property is not checked vacuously.
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>

namespace
{
    namespace r = nxx::roots;
    using nxx::stop_reason;

    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double inf = std::numeric_limits<double>::infinity();

    // Opposite signs or an exact zero, by comparisons.
    bool changes_sign(double a, double b) { return (a <= 0.0 && b >= 0.0) || (a >= 0.0 && b <= 0.0); }

    // abs + rel * m, rounded once per operation. Two statements, so that Clang does not contract them into an fma: the
    // library computes its thresholds without contraction, and the checks below are exact comparisons.
    double mixed_bound(double abs_part, double rel, double m)
    {
        const double rel_part = rel * m;
        return abs_part + rel_part;
    }

    // Brent's documented guarantee (roots/brent.hpp, DESIGN §6.8): width <= threshold(b) + 4 eps |b|.
    // brent's resolution_limit: the tolerance was below its resolution floor (2 eps |b| per side), so the enclosure is
    // at most 4 eps |b| wide (roots/brent.hpp).
    template<class E>
    bool brent_resolution(const E& e, double b)
    { return e.width() <= 4.0 * eps * std::abs(b); }

    // bisection's resolution_limit: the enclosure cannot be split, its ends are adjacent (DESIGN §9.3).
    struct adjacent_ends
    {
        template<class E>
        bool operator()(const E& e, double) const
        { return e.hi() == std::nextafter(e.lo(), inf); }
    };

    // floored_width{bits}: max(2^(1 - bits), 4 eps) (bits = 53 for the default), computed here independently.
    double floored_factor(int bits) { return (std::max)(std::ldexp(1.0, 1 - bits), 4.0 * eps); }

    // step_tol<Num, Den>: 2^-ceil(53 Num / Den) max(|x|, 1) for double.
    constexpr int step_exponent(int num, int den) { return (std::numeric_limits<double>::digits * num + den - 1) / den; }
    static_assert(step_exponent(7, 10) == 38);    // secant's default
    static_assert(step_exponent(3, 5) == 32);     // newton's default
    double step_bound(int num, int den, double x) { return std::ldexp(1.0, -step_exponent(num, den)) * (std::max)(std::abs(x), 1.0); }

    // f(x) = g(y), y = k (x - root) - s, increasing in x. `root` is a double and s = frac k gap with 0 < frac < 1 and
    // gap = nextafter(root, +inf) - root, so the true root lies strictly between the adjacent doubles root and
    // root + gap. Near it x - root is exact and k (x - root) never equals s, so for shapes 0 and 1 the sign of f(x) is
    // exactly the sign of x - x*, and no double is an exact zero: a bisection whose tolerance cannot be met ends at
    // the resolution limit, with the enclosure [root, root + gap]. Shape 2 (exp(y) - 1) is 0 wherever exp(y) rounds
    // to 1, so it also has exact zeros off the root.
    struct shaped_fn
    {
        int    shape = 0;      // 0: y (1 + c y^2), 1: tanh(y), 2: exp(y) - 1
        double root  = 0.0;    // the double just below the true root
        double k     = 1.0;    // 1 / the problem's scale
        double s     = 0.0;    // the offset that moves the true root off root
        double c     = 1.0;    // > 0: the curvature of shape 0

        double arg(double x) const { return k * (x - root) - s; }

        double operator()(double x) const
        {
            const double y = arg(x);
            switch (shape) {
                case 0:
                    return y * (1.0 + c * y * y);
                case 1:
                    return std::tanh(y);
                default:
                    return std::exp(y) - 1.0;
            }
        }

        double slope(double x) const
        {
            const double y = arg(x);
            switch (shape) {
                case 0:
                    return k * (1.0 + 3.0 * c * y * y);
                case 1: {
                    const double t = std::tanh(y);
                    return k * (1.0 - t * t);
                }
                default:
                    return k * std::exp(y);
            }
        }

        double above() const { return std::nextafter(root, inf); }    // the double just above the true root

        bool exact_sign() const { return shape != 2; }    // f(x) < 0 exactly for x <= root, f(x) > 0 for x >= above()

        // Whether [lo, hi] contains the true root (shapes 0 and 1 only).
        bool encloses_root(double lo, double hi) const { return lo <= root && above() <= hi; }
    };

    struct test_problem
    {
        shaped_fn fn;
        double    lo = 0.0;    // lo < root < hi
        double    hi = 0.0;
        double    x0 = 0.0;    // a guess near the root, for the open methods
    };

    class problem_source
    {
        std::mt19937 gen_;

        double uniform(double a, double b) { return std::uniform_real_distribution<double>(a, b)(gen_); }

    public:
        explicit problem_source(std::uint32_t seed) : gen_(seed) {}

        int integer(int a, int b) { return std::uniform_int_distribution<int>(a, b)(gen_); }

        double power_of_ten(int lo_exp, int hi_exp) { return std::pow(10.0, integer(lo_exp, hi_exp)); }

        // A root at a magnitude from 1e-8 to 1e8, either sign. Half the brackets stay on the root's side of 0 (scale
        // |root|); the others have scale max(|root|, 1), so they straddle 0 when |root| < 1.
        test_problem next()
        {
            test_problem p;
            p.fn.shape         = integer(0, 2);
            const double mag   = std::pow(10.0, uniform(-8.0, 8.0));
            const bool   neg   = integer(0, 1) == 1;
            const bool   local = integer(0, 1) == 1;
            const double base  = local ? mag : (std::max)(mag, 1.0);
            const double below = 0.9 * base * std::pow(10.0, uniform(-4.0, 0.0));
            const double above = 0.9 * base * std::pow(10.0, uniform(-4.0, 0.0));
            const double frac  = uniform(0.1, 0.9);
            p.fn.root          = neg ? -mag : mag;
            p.lo               = p.fn.root - below;
            p.hi               = p.fn.root + above;
            p.fn.k             = std::pow(10.0, uniform(-1.0, 1.0)) / (std::max)(below, above);
            p.fn.s             = frac * p.fn.k * (p.fn.above() - p.fn.root);
            p.fn.c             = std::pow(10.0, uniform(-2.0, 2.0));
            p.x0               = p.fn.root + uniform(-0.5, 0.5) / p.fn.k;
            return p;
        }
    };

    struct tally
    {
        int criterion        = 0;
        int exact_zero       = 0;
        int resolution_limit = 0;
        int failed           = 0;
    };

    // A bracketing solver's success: a sign-changing enclosure that contains the root, x at one of its ends with
    // fx == f(x), and the guarantee of the stop reason. width_ok(enclosure, x) is the criterion's guarantee.
    template<class R, class WidthOk, class ResolutionOk = adjacent_ends>
    void check_bracketing(const R& res, const test_problem& p, WidthOk width_ok, tally& t, ResolutionOk resolution_ok = {})
    {
        if (!res) {
            ++t.failed;
            return;
        }
        CHECK(res->fx == p.fn(res->x));
        if (!res->enclosure) {    // brent at b == c: a zero-width result, within every bound
            CHECK(res->uncertainty == 0.0);
            CHECK((res->how == stop_reason::criterion || res->how == stop_reason::exact_zero));
            if (res->how == stop_reason::criterion) ++t.criterion;
            if (res->how == stop_reason::exact_zero) {
                ++t.exact_zero;
                CHECK(res->fx == 0.0);
            }
            return;
        }
        const auto& e = *res->enclosure;
        CHECK(e.lo() < e.hi());
        CHECK(changes_sign(e.flo(), e.fhi()));
        CHECK(e.flo() == p.fn(e.lo()));
        CHECK(e.fhi() == p.fn(e.hi()));
        CHECK(res->uncertainty == e.width());
        CHECK((res->x == e.lo() || res->x == e.hi()));
        if (p.fn.exact_sign()) CHECK(p.fn.encloses_root(e.lo(), e.hi()));
        switch (res->how) {
            case stop_reason::criterion:
                ++t.criterion;
                CHECK(width_ok(e, res->x));
                break;
            case stop_reason::resolution_limit:
                ++t.resolution_limit;
                CHECK(resolution_ok(e, res->x));
                break;
            case stop_reason::exact_zero:
                ++t.exact_zero;
                CHECK(res->fx == 0.0);
                break;
            default:
                FAIL_CHECK("unexpected stop reason for a bracketing solver");
        }
    }

    // An open method's success: fx == f(x), no enclosure, and the guarantee of the stop reason.
    template<class R, class Ok>
    void check_open(const R& res, const test_problem& p, Ok ok, tally& t)
    {
        if (!res) {
            ++t.failed;
            return;
        }
        CHECK(res->fx == p.fn(res->x));
        CHECK_FALSE(res->enclosure.has_value());
        switch (res->how) {
            case stop_reason::criterion:
                ++t.criterion;
                CHECK(ok(*res));
                break;
            case stop_reason::exact_zero:
                ++t.exact_zero;
                CHECK(res->fx == 0.0);
                break;
            default:
                FAIL_CHECK("unexpected stop reason for an open method");
        }
    }

    auto derivative_of_problem(const shaped_fn& fn)
    {
        return [&fn](double x) { return fn.slope(x); };
    }

    // Whether s.with_stop(c) is well-formed (in a template, so that an ill-formed call makes it false).
    template<class S, class C>
    inline constexpr bool with_stop_ok_v = requires(const S& s, const C& c) { s.with_stop(c); };

    // Whether S{c} deduces and constructs a solver (class template argument deduction on the template name).
    template<class C>
    inline constexpr bool bisection_from_v = requires(const C& c) { r::bisection { c }; };
    template<class C>
    inline constexpr bool brent_from_v = requires(const C& c) { r::brent { c }; };
}    // namespace

TEST_SUITE("roots")
{
    TEST_CASE("soundness: bisection with width_tol on random problems")
    {
        problem_source src(20260928u);
        tally          t;
        for (int i = 0; i < 700; ++i) {
            const test_problem p = src.next();
            // abs in {0, 1e-12, ..., 1e-1}, rel in {0, 1e-15, ..., 1e-3}, not both 0.
            const int    mode = src.integer(0, 2);    // 0: abs only, 1: rel only, 2: both
            const double abs  = mode == 1 ? 0.0 : src.power_of_ten(-12, -1);
            const double rel  = mode == 0 ? 0.0 : src.power_of_ten(-15, -3);
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);
            CAPTURE(abs);
            CAPTURE(rel);
            const auto wt = nxx::rel_tolerance<double>::make(rel).and_then([&](auto rp) { return nxx::width_tol<double>::make(abs, rp); });
            if (!wt) {
                FAIL_CHECK("width_tol::make rejected a valid tolerance");
                continue;
            }
            const auto res = r::bisection { *wt }(p.fn, { p.lo, p.hi });
            check_bracketing(
                res,
                p,
                [&](const auto& e, double) {
                    return e.hi() - e.lo() <= mixed_bound(abs, rel, (std::min)(std::abs(e.lo()), std::abs(e.hi())));
                },
                t);
            if (res) CHECK(res->by == r::algos::bisection);
        }
        MESSAGE("criterion ", t.criterion, ", resolution_limit ", t.resolution_limit, ", exact_zero ", t.exact_zero, ", failed ", t.failed);
        CHECK(t.criterion >= 300);
        CHECK(t.resolution_limit >= 5);
    }

    TEST_CASE("soundness: bisection with floored_width on random problems")
    {
        problem_source src(1234567u);
        tally          t;
        for (int i = 0; i < 300; ++i) {
            const test_problem p = src.next();
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);
            const auto floored_ok = [](int bits) {
                return [bits](const auto& e, double) {
                    return e.hi() - e.lo() <= floored_factor(bits) * (std::max)(1.0, (std::min)(std::abs(e.lo()), std::abs(e.hi())));
                };
            };
            check_bracketing(r::bisection {}(p.fn, { p.lo, p.hi }), p, floored_ok(53), t);
            check_bracketing(r::bisection { nxx::floored_width {} }(p.fn, { p.lo, p.hi }), p, floored_ok(53), t);
            check_bracketing(r::bisection { nxx::floored_width { 20 } }(p.fn, { p.lo, p.hi }), p, floored_ok(20), t);
        }
        MESSAGE("criterion ", t.criterion, ", resolution_limit ", t.resolution_limit, ", exact_zero ", t.exact_zero, ", failed ", t.failed);
        CHECK(t.criterion >= 600);
    }

    TEST_CASE("soundness: brent meets the width criterion or ends at its resolution floor")
    {
        problem_source src(777u);
        tally          t;
        for (int i = 0; i < 400; ++i) {
            const test_problem p    = src.next();
            const int          mode = src.integer(0, 2);
            const double       abs  = mode == 1 ? 0.0 : src.power_of_ten(-12, -1);
            const double       rel  = mode == 0 ? 0.0 : src.power_of_ten(-15, -3);
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);
            CAPTURE(abs);
            CAPTURE(rel);
            const auto wt = nxx::rel_tolerance<double>::make(rel).and_then([&](auto rp) { return nxx::width_tol<double>::make(abs, rp); });
            if (!wt) {
                FAIL_CHECK("width_tol::make rejected a valid tolerance");
                continue;
            }
            const auto res_w = r::brent { *wt }(p.fn, { p.lo, p.hi });
            check_bracketing(
                res_w,
                p,
                [&](const auto& e, double) { return e.width() <= mixed_bound(abs, rel, (std::min)(std::abs(e.lo()), std::abs(e.hi()))); },
                t,
                [](const auto& e, double b) { return brent_resolution(e, b); });
            // The zero-width case has no enclosure: its uncertainty is 0, within every bound.
            if (res_w) CHECK(res_w->by == r::algos::brent);

            const auto res_f = r::brent {}(p.fn, { p.lo, p.hi });
            check_bracketing(
                res_f,
                p,
                [&](const auto& e, double) {
                    return e.width() <= floored_factor(53) * (std::max)(1.0, (std::min)(std::abs(e.lo()), std::abs(e.hi())));
                },
                t,
                [](const auto& e, double b) { return brent_resolution(e, b); });
        }
        MESSAGE("criterion ", t.criterion, ", resolution_limit ", t.resolution_limit, ", exact_zero ", t.exact_zero, ", failed ", t.failed);
        CHECK(t.criterion >= 600);
    }

    TEST_CASE("soundness: secant with x_tol and step_tol and f_tol")
    {
        problem_source src(4242u);
        tally          tx;
        tally          ts;
        tally          tf;
        for (int i = 0; i < 300; ++i) {
            const test_problem p    = src.next();
            const int          mode = src.integer(0, 2);
            const double       abs  = mode == 1 ? 0.0 : src.power_of_ten(-14, -2);
            const double       rel  = mode == 0 ? 0.0 : src.power_of_ten(-15, -4);
            const double       ftol = src.power_of_ten(-14, -2);
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.x0);
            CAPTURE(abs);
            CAPTURE(rel);
            CAPTURE(ftol);
            const auto xt = nxx::rel_tolerance<double>::make(rel).and_then([&](auto rp) { return nxx::x_tol<double>::make(abs, rp); });
            const auto ft = nxx::tolerance<double>::make(ftol);
            if (!xt || !ft) {
                FAIL_CHECK("make rejected a valid tolerance");
                continue;
            }
            check_open(
                r::secant { *xt }(p.fn, p.x0),
                p,
                [&](const auto& e) { return e.uncertainty <= mixed_bound(abs, rel, std::abs(e.x)); },
                tx);
            check_open(
                r::secant {}(p.fn, p.x0),
                p,
                [&](const auto& e) {
                    return e.uncertainty <= nxx::step_tol<7, 10>::threshold(e.x) && e.uncertainty <= step_bound(7, 10, e.x);
                },
                ts);
            check_open(r::secant { nxx::f_tol<double> { *ft } }(p.fn, p.x0), p, [&](const auto& e) { return std::abs(e.fx) <= ftol; }, tf);
        }
        MESSAGE("x_tol: criterion ", tx.criterion, ", exact_zero ", tx.exact_zero, ", failed ", tx.failed);
        MESSAGE("step_tol: criterion ", ts.criterion, ", exact_zero ", ts.exact_zero, ", failed ", ts.failed);
        MESSAGE("f_tol: criterion ", tf.criterion, ", exact_zero ", tf.exact_zero, ", failed ", tf.failed);
        CHECK(tx.criterion >= 100);
        CHECK(ts.criterion >= 100);
        CHECK(tf.criterion >= 100);
    }

    TEST_CASE("soundness: newton with x_tol and step_tol and f_tol")
    {
        problem_source src(99991u);
        tally          tx;
        tally          ts;
        tally          tf;
        for (int i = 0; i < 300; ++i) {
            const test_problem p    = src.next();
            const int          mode = src.integer(0, 2);
            const double       abs  = mode == 1 ? 0.0 : src.power_of_ten(-14, -2);
            const double       rel  = mode == 0 ? 0.0 : src.power_of_ten(-15, -4);
            const double       ftol = src.power_of_ten(-14, -2);
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.x0);
            CAPTURE(abs);
            CAPTURE(rel);
            CAPTURE(ftol);
            const auto xt = nxx::rel_tolerance<double>::make(rel).and_then([&](auto rp) { return nxx::x_tol<double>::make(abs, rp); });
            const auto ft = nxx::tolerance<double>::make(ftol);
            if (!xt || !ft) {
                FAIL_CHECK("make rejected a valid tolerance");
                continue;
            }
            const auto df = derivative_of_problem(p.fn);
            check_open(
                r::newton { *xt }.with_derivative(df)(p.fn, p.x0),
                p,
                [&](const auto& e) { return e.uncertainty <= mixed_bound(abs, rel, std::abs(e.x)); },
                tx);
            check_open(
                r::newton {}.with_derivative(df)(p.fn, p.x0),
                p,
                [&](const auto& e) {
                    return e.uncertainty <= nxx::step_tol<3, 5>::threshold(e.x) && e.uncertainty <= step_bound(3, 5, e.x);
                },
                ts);
            check_open(
                r::newton { nxx::f_tol<double> { *ft } }.with_derivative(df)(p.fn, p.x0),
                p,
                [&](const auto& e) { return std::abs(e.fx) <= ftol; },
                tf);
        }
        MESSAGE("x_tol: criterion ", tx.criterion, ", exact_zero ", tx.exact_zero, ", failed ", tx.failed);
        MESSAGE("step_tol: criterion ", ts.criterion, ", exact_zero ", ts.exact_zero, ", failed ", ts.failed);
        MESSAGE("f_tol: criterion ", tf.criterion, ", exact_zero ", tf.exact_zero, ", failed ", tf.failed);
        CHECK(tx.criterion >= 150);
        CHECK(ts.criterion >= 150);
        CHECK(tf.criterion >= 150);
    }

    TEST_CASE("soundness: f_tol on the bracketing solvers")
    {
        problem_source src(31337u);
        tally          tb;
        tally          tr;
        for (int i = 0; i < 400; ++i) {
            const test_problem p    = src.next();
            const double       ftol = src.power_of_ten(-14, -2);
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);
            CAPTURE(ftol);
            const auto ft = nxx::tolerance<double>::make(ftol);
            if (!ft) {
                FAIL_CHECK("tolerance::make rejected a valid tolerance");
                continue;
            }
            const nxx::f_tol<double> crit { *ft };
            check_bracketing(
                r::bisection {}.with_stop(crit)(p.fn, { p.lo, p.hi }),
                p,
                [&](const auto&, double x) { return std::abs(p.fn(x)) <= ftol; },
                tb);
            // brent keeps its width tolerance: a criterion stop meets the width guarantee or the f_tol.
            check_bracketing(
                r::brent {}.with_stop(crit)(p.fn, { p.lo, p.hi }),
                p,
                [&](const auto& e, double b) {
                    return std::abs(p.fn(b)) <= ftol ||
                           e.width() <= floored_factor(53) * (std::max)(1.0, (std::min)(std::abs(e.lo()), std::abs(e.hi())));
                },
                tr,
                [](const auto& e, double b) { return brent_resolution(e, b); });
        }
        MESSAGE("bisection: criterion ", tb.criterion, ", resolution_limit ", tb.resolution_limit, ", failed ", tb.failed);
        MESSAGE("brent: criterion ", tr.criterion, ", failed ", tr.failed);
        CHECK(tb.criterion >= 150);
        CHECK(tr.criterion >= 150);
    }

    TEST_CASE("soundness: bracketing solvers never converge while the width exceeds the tolerance near the resolution limit")
    {
        // Relative tolerances of a few ulps: the criterion is met only in the last halvings, where a criterion stop
        // and a resolution-limit stop compete.
        problem_source src(8675309u);
        tally          tb;
        tally          tr;
        for (int i = 0; i < 400; ++i) {
            const test_problem p   = src.next();
            const double       rel = std::ldexp(1.0, src.integer(-53, -48));    // 0.5 .. 16 eps
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);
            CAPTURE(rel);
            const auto wt = nxx::rel_tolerance<double>::make(rel).and_then([](auto rp) { return nxx::width_tol<double>::make(0.0, rp); });
            if (!wt) {
                FAIL_CHECK("width_tol::make rejected a valid tolerance");
                continue;
            }
            check_bracketing(
                r::bisection { *wt }(p.fn, { p.lo, p.hi }),
                p,
                [&](const auto& e, double) { return e.hi() - e.lo() <= rel * (std::min)(std::abs(e.lo()), std::abs(e.hi())); },
                tb);
            check_bracketing(
                r::brent { *wt }(p.fn, { p.lo, p.hi }),
                p,
                [&](const auto& e, double) { return e.width() <= rel * (std::min)(std::abs(e.lo()), std::abs(e.hi())); },
                tr,
                [](const auto& e, double b) { return brent_resolution(e, b); });
        }
        MESSAGE("bisection: criterion ",
                tb.criterion,
                ", resolution_limit ",
                tb.resolution_limit,
                ", exact_zero ",
                tb.exact_zero,
                ", failed ",
                tb.failed);
        MESSAGE("brent: criterion ",
                tr.criterion,
                ", resolution_limit ",
                tr.resolution_limit,
                ", exact_zero ",
                tr.exact_zero,
                ", failed ",
                tr.failed);
        CHECK(tb.criterion >= 50);
        CHECK(tb.resolution_limit >= 50);
        CHECK(tr.criterion >= 100);
        CHECK(tr.resolution_limit >= 50);    // tolerances below brent's resolution floor (4 eps |x|)
    }

    TEST_CASE("soundness: resolution_limit implies adjacent endpoints")
    {
        // width_tol{1e-300} is unreachable at |x| >= 1e-8, so every bisection ends at the resolution limit, or at an
        // exact zero of shape 2. For shapes 0 and 1 the final enclosure is exactly [root, nextafter(root, +inf)].
        problem_source src(55555u);
        int            limits = 0;
        for (int i = 0; i < 300; ++i) {
            const test_problem p = src.next();
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);
            const auto res = r::bisection { nxx::width_tol { 1e-300 } }(p.fn, { p.lo, p.hi });
            if (!res) {
                FAIL_CHECK("bisection with an unreachable tolerance failed");
                continue;
            }
            if (p.fn.exact_sign())
                CHECK(res->how == stop_reason::resolution_limit);
            else
                CHECK((res->how == stop_reason::resolution_limit || res->how == stop_reason::exact_zero));
            if (res->how != stop_reason::resolution_limit) continue;
            ++limits;
            if (res->enclosure) {
                CHECK(res->enclosure->hi() == std::nextafter(res->enclosure->lo(), inf));
                CHECK(res->uncertainty == res->enclosure->width());
                if (p.fn.exact_sign()) {
                    CHECK(res->enclosure->lo() == p.fn.root);
                    CHECK(res->enclosure->hi() == p.fn.above());
                }
            }
            else
                FAIL_CHECK("a resolution-limit result without an enclosure");
        }
        MESSAGE("resolution_limit ", limits);
        CHECK(limits >= 180);

        // Near x ~ 1, the case the design names: x - 2 / x on [1, 2], root sqrt(2), where no double is an exact zero
        // (2 / x rounds to the other neighbour of sqrt(2) at both of them).
        const auto sq  = [](double x) { return x - 2.0 / x; };
        const auto res = r::bisection { nxx::width_tol { 1e-300 } }(sq, { 1.0, 2.0 });
        if (res && res->enclosure) {
            CHECK(res->how == stop_reason::resolution_limit);
            CHECK(res->enclosure->hi() == std::nextafter(res->enclosure->lo(), inf));
            CHECK(res->enclosure->lo() < std::sqrt(2.0) + 1e-15);
            CHECK(res->enclosure->hi() > std::sqrt(2.0) - 1e-15);
        }
        else
            FAIL_CHECK("bisection on sqrt(2) failed or has no enclosure");
    }

    TEST_CASE("soundness: regression on exp minus 1.0001 with a 1e-3 width tolerance")
    {
        const auto   fn   = [](double x) { return std::exp(x) - 1.0001; };
        const double root = std::log1p(1.0001 - 1.0);    // 1.0001 - 1.0 is exact

        const auto bis = r::bisection { nxx::width_tol { 1e-3 } }(fn, { 0.0, 10.0 });
        if (bis && bis->enclosure) {
            CHECK(bis->how == stop_reason::criterion);
            CHECK(bis->enclosure->width() <= 1e-3);
            CHECK(bis->uncertainty <= 1e-3);
            CHECK(bis->enclosure->lo() <= root);
            CHECK(root <= bis->enclosure->hi());
            CHECK(std::abs(bis->x - root) <= 1e-3);
            CHECK(bis->fx == fn(bis->x));
        }
        else
            FAIL_CHECK("bisection failed or has no enclosure");

        const auto bre = r::brent { nxx::width_tol { 1e-3 } }(fn, { 0.0, 10.0 });
        if (bre) {
            CHECK(bre->how == stop_reason::criterion);
            const double bound = 1e-3;    // the width criterion itself (DESIGN §9.3)
            CHECK(bre->uncertainty <= bound);
            CHECK(std::abs(bre->x - root) <= bound);
            CHECK(bre->fx == fn(bre->x));
            if (bre->enclosure) {
                CHECK(bre->enclosure->width() <= bound);
                CHECK(bre->enclosure->lo() <= root);
                CHECK(root <= bre->enclosure->hi());
            }
        }
        else
            FAIL_CHECK("brent failed");
    }

    TEST_CASE("soundness: every intermediate bracket keeps a sign change and strictly shrinks")
    {
        problem_source src(2718281u);
        int            walked = 0;
        for (int i = 0; i < 150; ++i) {
            const test_problem p = src.next();
            CAPTURE(i);
            CAPTURE(p.fn.shape);
            CAPTURE(p.fn.root);
            CAPTURE(p.lo);
            CAPTURE(p.hi);

            // bisection: every state's sign_bracket.
            {
                const auto solver = r::bisection {};
                auto       prob   = solver.prepare(std::cref(p.fn), std::pair { p.lo, p.hi });
                if (!prob) {
                    FAIL_CHECK("bisection::prepare failed on a sign-changing bracket");
                    continue;
                }
                std::optional<r::sign_bracket<double>> prev;
                int                                    k = 0;
                for (const auto& st : nxx::steps_view { solver, *prob } | std::views::take(120)) {
                    CAPTURE(k);
                    if (!st) {
                        FAIL_CHECK("a bisection step failed");
                        break;
                    }
                    const r::sign_bracket<double>& b = st->b;
                    CHECK(b.lo() < b.hi());
                    CHECK(changes_sign(b.flo(), b.fhi()));
                    CHECK(b.flo() == p.fn(b.lo()));
                    CHECK(b.fhi() == p.fn(b.hi()));
                    if (p.fn.exact_sign()) CHECK(p.fn.encloses_root(b.lo(), b.hi()));
                    if (prev) {
                        CHECK(prev->lo() <= b.lo());
                        CHECK(b.hi() <= prev->hi());
                        CHECK(b.width() < prev->width());
                    }
                    prev = b;
                    ++k;
                }
                CHECK(k >= 2);
                CHECK(k < 120);    // ends at an intrinsic stop (exact zero or resolution limit)
                ++walked;
            }

            // brent: the enclosure of every state's estimate (none once b == c).
            {
                const auto solver = r::brent {};
                auto       prob   = solver.prepare(std::cref(p.fn), std::pair { p.lo, p.hi });
                if (!prob) {
                    FAIL_CHECK("brent::prepare failed on a sign-changing bracket");
                    continue;
                }
                std::optional<r::sign_bracket<double>> prev;
                int                                    k = 0;
                for (const auto& st : nxx::steps_view { solver, *prob } | std::views::take(120)) {
                    CAPTURE(k);
                    if (!st) {
                        FAIL_CHECK("a brent step failed");
                        break;
                    }
                    const auto e = solver.estimate(*st);
                    ++k;
                    if (!e.enclosure) {
                        CHECK(e.uncertainty == 0.0);
                        continue;
                    }
                    const r::sign_bracket<double>& b = *e.enclosure;
                    CHECK(b.lo() < b.hi());
                    CHECK(changes_sign(b.flo(), b.fhi()));
                    CHECK(b.flo() == p.fn(b.lo()));
                    CHECK(b.fhi() == p.fn(b.hi()));
                    if (p.fn.exact_sign()) CHECK(p.fn.encloses_root(b.lo(), b.hi()));
                    if (prev) {
                        CHECK(prev->lo() <= b.lo());
                        CHECK(b.hi() <= prev->hi());
                        CHECK(b.width() < prev->width());
                    }
                    prev = b;
                }
                CHECK(k >= 1);
                CHECK(k < 120);    // ends at an intrinsic stop (brent's tolerance or an exact zero)
            }
        }
        CHECK(walked == 150);
    }

    TEST_CASE("soundness: a width tolerance whose threshold overflows never accepts an infinite width")
    {
        // width_tol{max, nxx::rel_tolerance{0.1}} on [-max, max]: abs + rel min(|lo|, |hi|) = 1.1 max and the width 2 max were both
        // computed as inf, so brent reported criterion before its first step (inf <= inf), at x = max, 1.9 max from the root. The checks
        // compare halves, which cannot overflow.
        constexpr double big  = (std::numeric_limits<double>::max)();
        const double     root = -0.9 * big;
        const auto       fn   = [root](double x) {    // continuous and increasing; f(-max) = -1, f(max) = 0.5
            const double y = std::tanh(x / 1e306 - root / 1e306);
            return y > 0.0 ? 0.5 * y : y;
        };
        const auto wt = nxx::width_tol<double>::make(big, nxx::rel_tolerance { 0.1 });
        if (!wt) {
            FAIL_CHECK("width_tol::make rejected a valid tolerance");
            return;
        }
        const auto check_halves = [&](const auto& res, const char* name) {
            INFO(name);
            if (!res || !res->enclosure) {
                FAIL_CHECK("no success with an enclosure");
                return;
            }
            const auto&  e    = *res->enclosure;
            const double half = big / 2.0 + 0.1 * (std::min)(std::abs(e.lo()), std::abs(e.hi())) / 2.0;    // (abs + rel m) / 2
            CHECK(res->how == stop_reason::criterion);
            CHECK(res->used.iterations >= 1u);
            CHECK(std::isfinite(res->uncertainty));
            CHECK(e.hi() / 2.0 - e.lo() / 2.0 <= half);
            CHECK(std::abs(res->x / 2.0 - root / 2.0) <= half);
            CHECK(e.lo() <= root);
            CHECK(root <= e.hi());
            CHECK(res->fx == fn(res->x));
        };
        check_halves(r::brent { *wt }(fn, { -big, big }), "brent");
        check_halves(r::bisection { *wt }(fn, { -big, big }), "bisection");    // tests only after a halving

        // A double tolerance above FLT_MAX on a float problem: float(5e38) is inf. The bound is checked in double, exactly.
        constexpr float fbig = (std::numeric_limits<float>::max)();
        const auto      ffn  = [](float x) { return std::tanh(x / 1e36f + 300.0f); };    // root near -3e38; f(+-max) = +-1
        const auto      fres = r::brent { nxx::width_tol { 5e38 } }(ffn, { -fbig, fbig });
        if (fres && fres->enclosure) {
            CHECK(fres->how == stop_reason::criterion);
            CHECK(static_cast<double>(fres->enclosure->hi()) - static_cast<double>(fres->enclosure->lo()) <= 5e38);
            CHECK(std::abs(static_cast<double>(fres->x) + 3e38) <= 5e38);
        }
        else
            FAIL_CHECK("brent on a float problem with width_tol{5e38} failed or has no enclosure");
    }

    TEST_CASE("soundness: x_tol on a bracketing solver does not compile")
    {
        using xt = nxx::x_tol<double>;
        using wt = nxx::width_tol<double>;
        using ft = nxx::f_tol<double>;

        // The constructors: deleted with a reason for x_tol, whatever the options type.
        static_assert(!std::is_constructible_v<r::bisection<>, xt>);
        static_assert(!std::is_constructible_v<r::bisection<nxx::options<xt>>, xt>);
        static_assert(!std::is_constructible_v<r::brent<>, xt>);
        static_assert(!std::is_constructible_v<r::brent<xt>, xt>);
        static_assert(!bisection_from_v<xt>);
        static_assert(!brent_from_v<xt>);
        // with_stop: the deleted sibling.
        static_assert(!with_stop_ok_v<r::bisection<>, xt>);
        static_assert(!with_stop_ok_v<r::brent<>, xt>);
        static_assert(!with_stop_ok_v<r::bisection<>, nxx::step_tol<7, 10>>);
        static_assert(!with_stop_ok_v<r::bisection<>, decltype(std::declval<xt>() || std::declval<ft>())>);

        // Controls: the width criteria and f_tol are accepted.
        static_assert(std::is_constructible_v<r::bisection<nxx::options<wt>>, wt>);
        static_assert(std::is_constructible_v<r::bisection<>, nxx::floored_width>);
        static_assert(std::is_constructible_v<r::brent<wt>, wt>);
        static_assert(bisection_from_v<wt>);
        static_assert(brent_from_v<wt>);
        static_assert(bisection_from_v<ft>);
        static_assert(with_stop_ok_v<r::bisection<>, wt>);
        static_assert(with_stop_ok_v<r::bisection<>, ft>);
        static_assert(with_stop_ok_v<r::brent<>, ft>);
        static_assert(with_stop_ok_v<r::bisection<>, decltype(std::declval<wt>() || std::declval<ft>())>);

        // And the mirror image: width criteria on the open methods do not compile; x_tol does.
        static_assert(!std::is_constructible_v<r::secant<nxx::options<wt>>, wt>);
        static_assert(!std::is_constructible_v<r::newton<nxx::options<wt>>, wt>);
        static_assert(!with_stop_ok_v<r::secant<>, wt>);
        static_assert(!with_stop_ok_v<r::newton<>, nxx::floored_width>);
        static_assert(std::is_constructible_v<r::secant<nxx::options<xt>>, xt>);
        static_assert(with_stop_ok_v<r::newton<>, xt>);
        CHECK(true);
    }
}
