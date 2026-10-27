// search(fn, x0) in nxx::roots (Phase 5, docs/roots/PHASE5_API.md): an outward search for a bracket.
//
// Evaluation counts and bracket ends are asserted where no libm is involved (+ - * / and sqrt in double), against an
// independent Python implementation of the specification (.scratch/roots5/proto5d.py), assuming no FMA contraction
// in fn, as in Phases 2 to 4. The ends are written as hexadecimal literals, exactly as the prototype prints them.

#include "roots_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <expected>
#include <limits>
#include <utility>

using namespace roots_check;

namespace {

    constexpr auto sq = [](double x) { return x * x - 2.0; };
    constexpr auto wal = [](double x) { return (x * x - 2.0) * x - 5.0; };
    constexpr auto reciprocal = [](double x) { return 1.0 / x - 1.0; };
    constexpr auto vdw = [](double v) { return 8.314 * 280.0 / (v - 4.27e-5) - 0.364 / (v * v) - 4e6; };
    // van der Waals with NaN outside its domain, v > b: the documented idiom.
    constexpr auto vdw_domain = [](double v) { return v > 4.27e-5 ? vdw(v) : dnan; };

    // search(fn, x0) with its evaluations counted.
    template<typename F>
    std::pair<std::expected<rt::Bracket<double>, rt::error>, int> searched(F fn, double x0)
    {
        int n = 0;
        auto const b = rt::search(Counted<F>{ fn, &n }, x0);
        return { b, n };
    }

    // fn's signs at the bracket's ends differ, or fn is 0 at one.
    template<typename F>
    bool sign_change(F const &fn, rt::Bracket<double> const &b)
    {
        double const lo = fn(b.lower);
        double const hi = fn(b.upper);
        return b.lower <= b.upper && (lo == 0.0 || hi == 0.0 || (lo < 0.0) != (hi < 0.0));
    }

} // namespace

// Constant evaluation in each type, and a domain edge (fn returns NaN explicitly: sqrt or log of a negative number is
// not a constant expression).
static_assert([] {
    auto const b = rt::search(sq, 1.0);
    return b && b->lower < sqrt2 && sqrt2 < b->upper;
}());
static_assert([] {
    auto const b = rt::search([](float x) { return x * x - 2.0f; }, 1.0f);
    return b && b->lower < 1.41421356f && 1.41421356f < b->upper;
}());
static_assert([] {
    auto const b = rt::search([](long double x) { return x * x - 2.0L; }, 1.0L);
    return b && b->lower < 1.41421356237309504880L && 1.41421356237309504880L < b->upper;
}());
static_assert([] {
    auto const b = rt::search([](double x) { return x > 0.0 ? x * x - 1e-6 : dnan; }, 1.0);
    return b && b->lower < 1e-3 && 1e-3 < b->upper;
}());
static_assert(rt::search(sq, dnan) == std::unexpected(rt::error::non_finite));
static_assert(rt::search([](double) { return 1.0; }, 1.0) == std::unexpected(rt::error::no_sign_change));

TEST_CASE("roots:search: evaluation counts and brackets in double")
{
    auto const check = [](auto const &fn, double x0, int count, double lower, double upper) {
        auto const [b, n] = searched(fn, x0);
        REQUIRE(b.has_value());
        CHECK(n == count);
        CHECK(b->lower == lower);
        CHECK(b->upper == upper);
    };
    check(sq, 1.0, 19, 0x1.577ae147ae148p+0, 0x1.833851eb851ecp+0);
    check(sq, 0.0, 24, -0x1.bade147ae147cp+0, -0x1.273eb851eb853p+0); // the step below comes first
    check(wal, 0.0, 27, 0x1.bade147ae147cp+0, 0x1.4c268f5c28f5dp+1);
    check(reciprocal, 3.0, 20, 0x1.63051eb851eb4p-1, 0x1.76570a3d70a3cp+0); // growth 2 would step over root and pole
    check([](double x) { return x - 1e6; }, 1.0, 91, 0x1.6c705dd02fdeap+19, 0x1.11543e5c23e70p+20);
    check([](double x) { return x - 1e-300; }, 1.0, 22, -0x1.39f5c28f5c298p-3, 0x1.d95c28f5c28f0p-3);
    // The first step |x0| / 50 underflows to 0 here; it is denorm_min instead.
    check([](double x) { return x + 1e-300; }, dmin, 266, -0x1.724e945583547p-997, -0x1.edbe1b1caf1b4p-998);
    check([](double x) { return x - 1.0; }, dmax / 2, 22, -0x1.39f5c28f5c28cp+1020, 0x1.d95c28f5c28f8p+1020);
    // A domain edge at 0 approached by bisection.
    check([](double x) { return std::sqrt(x) - 1e-3; }, 1.0, 60, 0x1.3f3d70a0181f8p-22, 0x1.14a3d7094db60p-20);
    check(vdw, 1e-3, 20, 0x1.e4b87bdcf02fcp-13, 0x1.ff1950331e3a4p-12);
    // Finite below its co-volume b, van der Waals has a pole at b, which the search brackets.
    check(vdw, 1e-4, 20, 0x1.83c6c97d8cf3cp-16, 0x1.98e10cf5b1c88p-15);
    check(vdw, 5e-5, 12, 0x1.63baba7b9170ep-15, 0x1.78f68be2f7b18p-15);
    // With NaN at v <= b, the search approaches b and finds the root.
    check(vdw_domain, 1e-4, 29, 0x1.78f49491f2dc3p-12, 0x1.00808c8259e20p-11);
    check(vdw_domain, 5e-5, 33, 0x1.6689efd86a8fdp-12, 0x1.ff9804d983948p-12);

    // No sign change: two roots within one step, none, an even root, a root beyond the cap. 2 cap + 1 evaluations.
    for (auto const &[x0, fn] : { std::pair{ 100.0, +[](double x) { return x * x - 2.0; } },
             std::pair{ 1.0, +[](double x) { return x * x + 1.0; } },
             std::pair{ 0.0, +[](double x) { return (x - 1.0) * (x - 1.0); } },
             std::pair{ 1e-300, +[](double x) { return x - 1.0; } } }) {
        auto const [b, n] = searched(fn, x0);
        CHECK(b == std::unexpected(rt::error::no_sign_change));
        CHECK(n == 425);
    }
}

TEST_CASE("roots:search: libm, the composition and the edges")
{
    // A domain edge through libm's -inf and NaN: ln(x / 1e-12) from 1.
    auto const ln = [](double x) { return std::log(x / 1e-12); };
    auto const b = rt::search(ln, 1.0);
    REQUIRE(b.has_value());
    CHECK(b->lower < 1e-12);
    CHECK(1e-12 < b->upper);
    CHECK(sign_change(ln, *b));

    // The composition with solve: the root (Phase 1's accuracy term for x^2 - 2), and an error passing through.
    auto const solved = [](auto const &fn, double x0) {
        return rt::search(fn, x0).and_then([&](rt::Bracket<double> br) { return rt::solve<rt::Brent>(fn, br); });
    };
    CHECK(brackets(solved(sq, 1.0), sqrt2, 8e-17));
    CHECK(solved([](double x) { return x * x + 1.0; }, 1.0) == std::unexpected(rt::error::no_sign_change));

    // fn exactly 0 at x0: {x0, x0}, which solve returns at once.
    auto const at_root = rt::search([](double x) { return x - 2.0; }, 2.0);
    CHECK(at_root == rt::Bracket<double>{ 2.0, 2.0 });
    CHECK(rt::solve<rt::Brent>([](double x) { return x - 2.0; }, *at_root) == rt::Root<double>{ 2.0, default_tol(2.0) });

    // x0 or fn(x0) not finite.
    CHECK(rt::search(sq, dinf) == std::unexpected(rt::error::non_finite));
    CHECK(rt::search([](double x) { return std::log(x); }, -1.0) == std::unexpected(rt::error::non_finite));

    // x0 at +-max: no overflow. No step exceeds max / (3/2), so from max the search reaches down to about max / 3: a
    // root at max / 2 is found, one at 1 is not.
    auto const half = [](double x) { return x - dmax / 2.0; };
    auto const neg_half = [](double x) { return x + dmax / 2.0; };
    auto const from_max = rt::search(half, dmax);
    auto const from_min = rt::search(neg_half, -dmax);
    REQUIRE(from_max.has_value());
    REQUIRE(from_min.has_value());
    CHECK(sign_change(half, *from_max));
    CHECK(sign_change(neg_half, *from_min));
    CHECK(rt::search([](double x) { return x - 1.0; }, dmax) == std::unexpected(rt::error::no_sign_change));

    // A tuple fn gives the value fn's bracket.
    auto const pair_fn = [](double x) { return std::pair{ (x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0 }; };
    CHECK(rt::search(pair_fn, 0.0) == rt::search(wal, 0.0));
}

TEST_CASE("roots:search: multiprecision types")
{
    auto const check = [](auto zero) {
        using T = decltype(zero);
        T const s2 = from_string<T>(sqrt2_70);
        auto const b = rt::search([](T const &x) { return T(x * x - 2); }, T(1));
        REQUIRE(b.has_value());
        CHECK(b->lower < s2);
        CHECK(s2 < b->upper);
    };
    check(bin50_et{});
    check(dec50{});

    // The cap in cpp_bin_float_50: 4 x 168 = 672 steps, 2 x 672 + 1 evaluations.
    int n = 0;
    auto const none = rt::search(
        [&n](bin50 const &x) {
            ++n;
            return bin50(x * x + 1);
        },
        bin50(1));
    CHECK(none == std::unexpected(rt::error::no_sign_change));
    CHECK(n == 1345);
}
