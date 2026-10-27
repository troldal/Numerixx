// The bracketing methods of nxx::roots, Bisection and Brent (Phase 1, docs/roots/PHASE1_API.md).
//
// References are closed forms and published constants, at 50 digits, cross-checked with mpmath (never Numerixx). A
// result must lie within abs_error of the reference plus a conditioning term c: the shift of fn's sign change from the
// exact root, fn's error budget near the root divided by |f'(root)|, derived per case (an fn of + - * /: its own
// rounding; a libm fn: a budget of 2 ulps). Evaluation counts are asserted only where no libm is involved (+ - * / in
// double), against an independent Python implementation of the specification (.scratch/roots1/proto3.py and
// counts.py; its values are in the comments); for Bisection they also follow from the derivation in the spec.

#include "roots_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>

using namespace roots_check;

namespace {

    constexpr auto square_minus_2 = [](double x) { return x * x - 2.0; };
    constexpr auto square_minus_4 = [](double x) { return x * x - 4.0; };
    constexpr auto minus_0375 = [](double x) { return x - 0.375; };
    constexpr auto cube = [](double x) { return x * x * x; };
    constexpr auto wallis_fn = [](double x) { return (x * x - 2.0) * x - 5.0; };
    constexpr auto minus_1 = [](double x) { return x - 1.0; };
    constexpr auto nan_below_0 = [](double x) { return x < 0.0 ? dnan : x - 1.0; };

} // namespace

// Exact cases in constant expressions: dyadic roots that a method reaches exactly.
static_assert(rt::solve<rt::Bisection>(minus_0375, rt::Bracket{ 0.0, 1.0 })->value == 0.375); // the third midpoint
static_assert(rt::solve<rt::Bisection>(square_minus_4, rt::Bracket{ 0.0, 4.0 })->value == 2.0); // the first midpoint
static_assert(rt::solve<rt::Brent>(square_minus_4, rt::Bracket{ 0.0, 4.0 })->value == 2.0);
// Zero first: a root at 0 in one evaluation more, reported with the tolerance at 0, denorm_min.
static_assert(rt::solve<rt::Brent>(cube, rt::Bracket{ -1.0, 2.0 }) == rt::Root<double>{ 0.0, dmin });
static_assert(rt::solve<rt::Bisection>(cube, rt::Bracket{ -1.0, 2.0 }) == rt::Root<double>{ 0.0, dmin });
// No overflow in constant evaluation, from the widest bracket (zero first gives [0, max], Brent's guard bisects).
static_assert(rt::solve<rt::Brent>(minus_1, rt::Bracket{ -dmax, dmax })->value == 1.0);
static_assert(nxx::math::abs(rt::solve<rt::Brent>(wallis_fn, rt::Bracket{ 2.0, 3.0 })->value - wallis) <= 1e-15);
static_assert(
    nxx::math::abs(rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 })->value - sqrt2) <= 1e-15);
static_assert(
    nxx::math::abs(
        rt::solve<rt::Brent>([](float x) { return x * x - 2.0f; }, rt::Bracket{ 1.0f, 2.0f })->value - 1.41421356f)
    <= 1e-6f);
static_assert(
    nxx::math::abs(rt::solve<rt::Brent>([](long double x) { return x * x - 2.0L; }, rt::Bracket{ 1.0L, 2.0L })->value
                   - 1.41421356237309504880L)
    <= 1e-15L);

// The independent review's overflow cases (Phase 1, 02fffe5) evaluate in constant expressions, with abs_error
// saturated at max: a tolerance sum at the edge of overflow, and a distance between ends of opposite signs whose
// overflow-free test had rounded (zero first skipped, since fn is NaN at 0).
static_assert(rt::solve<rt::Brent>([](double) { return 0.0; },
                  rt::Bracket{ 0x1.8p1023, 0x1.8p1023 },
                  { .abs_tol = 0x1.ffffffffffffep1021, .rel_tol = 1.0 })
              == rt::Root<double>{ 0x1.8p1023, dmax });
static_assert(rt::solve<rt::Brent>([](double x) { return x == 0.0 ? dnan : (x < 0.0 ? -1.0 : 2.0); },
                  rt::Bracket{ -0x1.0000000000003p1022, 0x1.7fffffffffffep1023 },
                  { .abs_tol = dmax / 4.0, .rel_tol = 1.0 })
                  ->abs_error
              == dmax);

// Failures are values in constant expressions too.
static_assert(
    rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 2.0, 3.0 }) == std::unexpected(rt::error::no_sign_change));
static_assert(
    rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ dnan, 2.0 }) == std::unexpected(rt::error::non_finite));
static_assert(rt::solve<rt::Brent>(nan_below_0, rt::Bracket{ -1.0, 4.0 }) == std::unexpected(rt::error::non_finite));
static_assert(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, { .abs_tol = -1.0 })
              == std::unexpected(rt::error::non_finite));

// The methods accept only their start, a Bracket of a real type.
template<typename M, typename S>
concept solvable = requires(S s) { rt::solve<M>(square_minus_2, s); };
static_assert(solvable<rt::Brent, rt::Bracket<double>> && solvable<rt::Bisection, rt::Bracket<double>>);
static_assert(!solvable<rt::Brent, double>); // a guess is not a bracket
static_assert(!solvable<int, rt::Bracket<double>>); // not a method
template<typename T>
concept bracket_of = requires { typename rt::Bracket<T>; };
static_assert(!bracket_of<int>); // integers are not real

TEST_CASE("roots: accuracy in double")
{
    // x^2 - 2: fn errs by at most half an ulp of 2 (x * x rounds; - 2 is exact near 2), so its sign change is within
    // 2.2e-16 / |f'| = 2.2e-16 / 2.83 = 7.9e-17 of sqrt(2).
    CHECK(brackets(rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }), sqrt2, 8e-17));
    CHECK(brackets(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }), sqrt2, 8e-17));

    // x^3 - 2x - 5 as (x^2 - 2) x - 5: about 19 half-ulp roundings of values up to 5, 9.6 eps, over f' = 11.2: 2e-16.
    CHECK(brackets(rt::solve<rt::Bisection>(wallis_fn, rt::Bracket{ 2.0, 3.0 }), wallis, 2e-16));
    CHECK(brackets(rt::solve<rt::Brent>(wallis_fn, rt::Bracket{ 2.0, 3.0 }), wallis, 2e-16));

    // cos x - x: 2 ulps of cos near 0.74 (2.2e-16; the difference is exact) over |f'| = 1 + sin(0.74) = 1.67: 1.4e-16.
    auto const dottie_fn = [](double x) { return std::cos(x) - x; };
    CHECK(brackets(rt::solve<rt::Bisection>(dottie_fn, rt::Bracket{ 0.0, 1.0 }), dottie, 1.4e-16));
    CHECK(brackets(rt::solve<rt::Brent>(dottie_fn, rt::Bracket{ 0.0, 1.0 }), dottie, 1.4e-16));

    // e^x - 2: 2 ulps of exp near 2 (8.9e-16) over f' = 2: 4.5e-16. Both methods meet an exact zero of fn here.
    auto const ln2_fn = [](double x) { return std::exp(x) - 2.0; };
    CHECK(brackets(rt::solve<rt::Bisection>(ln2_fn, rt::Bracket{ 0.0, 1.0 }), ln2, 4.5e-16));
    CHECK(brackets(rt::solve<rt::Brent>(ln2_fn, rt::Bracket{ 0.0, 1.0 }), ln2, 4.5e-16));

    // x e^x - 1: 2 ulps of exp near 1.76 times 0.57 (2.5e-16) and the product's half ulp (1.1e-16), over
    // f' = e^r (1 + r) = 2.76: 1.3e-16.
    auto const omega_fn = [](double x) { return x * std::exp(x) - 1.0; };
    CHECK(brackets(rt::solve<rt::Bisection>(omega_fn, rt::Bracket{ 0.0, 1.0 }), omega, 1.3e-16));
    CHECK(brackets(rt::solve<rt::Brent>(omega_fn, rt::Bracket{ 0.0, 1.0 }), omega, 1.3e-16));

    // A trace-level root keeps its relative accuracy with the default criteria: ln(x / 1e-12) changes sign within
    // about eps relative of fl(1e-12) (the division's and log's roundings), 2.2e-28.
    auto const trace_fn = [](double x) { return std::log(x / 1e-12); };
    CHECK(brackets(rt::solve<rt::Bisection>(trace_fn, rt::Bracket{ 1e-20, 1.0 }), 1e-12, 2.5e-28));
    CHECK(brackets(rt::solve<rt::Brent>(trace_fn, rt::Bracket{ 1e-20, 1.0 }), 1e-12, 2.5e-28));
}

TEST_CASE("roots: evaluation counts in double")
{
    // Bisection: 2 evaluations, then one per halving until half the bracket is at most tol(m): for [1, 2],
    // 2^-(k+1) <= 2 eps sqrt(2) at k = 50, so 52.
    auto const x20 = [](double x) {
        double const x2 = x * x;
        double const x4 = x2 * x2;
        double const x8 = x4 * x4;
        return x8 * x8 * x4 - 1.0;
    };
    auto const cubed = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };
    auto const tiny = [](double x) { return x - 1e-300; };

    CHECK(counted<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }).second == 52);
    CHECK(counted<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }).second == 9);
    CHECK(counted<rt::Bisection>(wallis_fn, rt::Bracket{ 2.0, 3.0 }).second == 51);
    CHECK(counted<rt::Brent>(wallis_fn, rt::Bracket{ 2.0, 3.0 }).second == 8);
    CHECK(counted<rt::Bisection>(x20, rt::Bracket{ 0.0, 5.0 }).second == 55);
    CHECK(counted<rt::Brent>(x20, rt::Bracket{ 0.0, 5.0 }).second == 19);
    // A triple root makes Brent's interpolation linear: slower than Bisection, as documented.
    CHECK(counted<rt::Bisection>(cubed, rt::Bracket{ 0.0, 3.0 }).second == 54);
    CHECK(counted<rt::Brent>(cubed, rt::Bracket{ 0.0, 3.0 }).second == 156);
    // A root at 1e-300 from [0, 1]: Bisection halves down to it, Brent's secant lands on it.
    CHECK(counted<rt::Bisection>(tiny, rt::Bracket{ 0.0, 1.0 }).second == 1049);
    CHECK(counted<rt::Brent>(tiny, rt::Bracket{ 0.0, 1.0 }).second == 3);
    // Brent's guard bisects while the bracket is wider than max / 16, then interpolates.
    CHECK(counted<rt::Bisection>(minus_1, rt::Bracket{ -dmax, dmax }).second == 1077);
    CHECK(counted<rt::Brent>(minus_1, rt::Bracket{ -dmax, dmax }).second == 8);
    CHECK(counted<rt::Brent>(minus_1, rt::Bracket{ 0.0, dmax }).second == 7);
    // Brent's worst case is (k + 1)^2 - 2 evaluations for k bisection steps (Brent 1973, ch. 4, eq. (3.4)).
    CHECK(counted<rt::Brent>(cubed, rt::Bracket{ 0.0, 3.0 }).second <= 53 * 53 - 2);
}

TEST_CASE("roots: the bracket")
{
    // Either order, the same root.
    CHECK(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 2.0, 1.0 })
          == rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }));
    CHECK(rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 2.0, 1.0 })
          == rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }));

    // A root at an end, and equal ends.
    CHECK(rt::solve<rt::Brent>(square_minus_4, rt::Bracket{ 2.0, 3.0 }) == rt::Root<double>{ 2.0, default_tol(2.0) });
    CHECK(
        rt::solve<rt::Brent>(minus_0375, rt::Bracket{ 0.375, 0.375 }) == rt::Root<double>{ 0.375, default_tol(0.375) });
    CHECK(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 1.0 }) == std::unexpected(rt::error::no_sign_change));

    // The same sign at both ends, even around two roots.
    CHECK(rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ -2.0, 2.0 })
          == std::unexpected(rt::error::no_sign_change));

    // Non-finite ends and end values.
    CHECK(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, dinf }) == std::unexpected(rt::error::non_finite));
    CHECK(rt::solve<rt::Brent>([](double x) { return 1.0 / x - 1.0; }, rt::Bracket{ 0.0, 2.0 })
          == std::unexpected(rt::error::non_finite)); // fn(0) = dinf
    CHECK(rt::solve<rt::Bisection>(nan_below_0, rt::Bracket{ -1.0, 4.0 }) == std::unexpected(rt::error::non_finite));
}

TEST_CASE("roots: fn not finite inside the bracket")
{
    // A NaN where the method evaluates: non_finite. Bisection's first midpoint is 0.5.
    auto const hole = [](double x) { return nxx::math::abs(x - 0.5) < 0.1 ? dnan : x - 0.5; };
    CHECK(rt::solve<rt::Bisection>(hole, rt::Bracket{ 0.0, 1.0 }) == std::unexpected(rt::error::non_finite));
    CHECK(rt::solve<rt::Brent>(hole, rt::Bracket{ 0.0, 1.0 }) == std::unexpected(rt::error::non_finite));

    // Not finite at 0 only, a removable singularity: zero first is skipped, and the root found. sin(x)/x - 0.5 has
    // its root at 1.8955 (mpmath); 2 ulps of sin over |f'| = 0.48 at the root give c = 5e-16. x / (e^x - 1) - 0.5
    // at 1.2564 (mpmath), with 2 ulps of expm1 and the quotient's rounding over |f'| = 0.32: c = 1.4e-15.
    auto const sinc = [](double x) { return std::sin(x) / x - 0.5; };
    auto const planck = [](double x) { return x / std::expm1(x) - 0.5; };
    CHECK(brackets(rt::solve<rt::Bisection>(sinc, rt::Bracket{ -1.0, 3.0 }), 1.8954942670339809471, 5e-16));
    CHECK(brackets(rt::solve<rt::Brent>(sinc, rt::Bracket{ -1.0, 3.0 }), 1.8954942670339809471, 5e-16));
    CHECK(brackets(rt::solve<rt::Bisection>(planck, rt::Bracket{ -1.0, 3.0 }), 1.2564312086261696770, 1.4e-15));
    CHECK(brackets(rt::solve<rt::Brent>(planck, rt::Bracket{ -1.0, 3.0 }), 1.2564312086261696770, 1.4e-15));
}

TEST_CASE("roots: criteria")
{
    // Negative and non-finite tolerances are errors.
    CHECK(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = -1e-8 })
          == std::unexpected(rt::error::non_finite));
    CHECK(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = dnan })
          == std::unexpected(rt::error::non_finite));
    CHECK(rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, { .abs_tol = dinf })
          == std::unexpected(rt::error::non_finite));

    // Below the floor they are raised to it: the same result as the defaults.
    rt::Criteria<double> const zero{ .abs_tol = 0.0, .rel_tol = 0.0 };
    CHECK(rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, zero)
          == rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }));
    CHECK(rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, zero)
          == rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }));
    CHECK(zero != rt::Criteria<double>{});

    // A looser tolerance: fewer evaluations, and the result within the bracket it asked for.
    auto const [loose, n_loose] = counted<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = 1e-6 });
    CHECK(loose.has_value());
    CHECK(n_loose < 9);
    CHECK(nxx::math::abs(loose->value - sqrt2) <= loose->abs_error + 8e-17);
    CHECK(loose->abs_error <= 2.0 * (dmin + 1e-6 * loose->value));

    // Huge tolerances are capped, and stop at once without overflow: one step, inside the bracket.
    for (rt::Criteria<double> const huge :
        { rt::Criteria<double>{ .rel_tol = 1e300 }, rt::Criteria<double>{ .abs_tol = dmax } }) {
        auto const b = rt::solve<rt::Bisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, huge);
        auto const r = rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }, huge);
        CHECK((b && b->value >= 1.0 && b->value <= 2.0 && b->abs_error <= 1.0));
        CHECK((r && r->value >= 1.0 && r->value <= 2.0 && r->abs_error <= 1.0));
    }

    // Loose tolerances stop with ends far apart, where their difference rounds: abs_error must still bound the distance
    // (the independent review's cases, 02fffe5, where it fell short by 0.5 and 0.25). The sign changes are at -0.5,
    // with fn NaN at 0 so that zero first is skipped, and at 1.75; the distances are checked exactly, in bin50.
    auto const at_minus_half = [](double x) { return x == 0.0 ? dnan : (x < -0.5 ? -2.0 : 1.0); };
    auto const at_175 = [](double x) { return x < 1.75 ? -1.0 : 1.0; };
    rt::Criteria<double> const loosest{ .rel_tol = 1.0 };
    CHECK(brackets_exactly(rt::solve<rt::Bisection>(at_minus_half, rt::Bracket{ -1.0, 0x1p54 }, loosest), bin50(-0.5)));
    CHECK(brackets_exactly(rt::solve<rt::Brent>(at_minus_half, rt::Bracket{ -1.0, 0x1p54 }, loosest), bin50(-0.5)));
    CHECK(brackets_exactly(rt::solve<rt::Bisection>(at_175, rt::Bracket{ 1.5, 0x1.8p54 }, loosest), bin50(1.75)));
    CHECK(brackets_exactly(rt::solve<rt::Brent>(at_175, rt::Bracket{ 1.5, 0x1.8p54 }, loosest), bin50(1.75)));

    // abs_tol near 0: a root at 1e-300 to an absolute 1e-12 in a few halvings rather than 1047.
    auto const [near0, n_near0] =
        counted<rt::Bisection>([](double x) { return x - 1e-300; }, rt::Bracket{ 0.0, 1.0 }, { .abs_tol = 1e-12 });
    CHECK((near0 && near0->value <= 2e-12));
    CHECK(n_near0 == 2 + 39); // 2^-(k+1) <= 1e-12 at k = 39
}

TEST_CASE("roots: edges of T")
{
    // Among the subnormals, where half the bracket rounds: 2x - 5 denorm_min on [0, 3 denorm_min], root 2.5
    // denorm_min. Every operation is exact, so the sign change is the root: c = 0.
    auto const sub = [](double x) { return 2.0 * x - 5.0 * dmin; };
    // 2.5 denorm_min is not a double (it rounds to 2 denorm_min), so the distance is checked exactly, in bin50.
    bin50 const sub_root = bin50(bin50(dmin) * bin50(5) / bin50(2));
    CHECK(brackets_exactly(rt::solve<rt::Bisection>(sub, rt::Bracket{ 0.0, 3.0 * dmin }), sub_root));
    CHECK(brackets_exactly(rt::solve<rt::Brent>(sub, rt::Bracket{ 0.0, 3.0 * dmin }), sub_root));

    // The widest brackets: [-max, max] (zero first gives [0, max]) and [0, max]; x - 1 is exact near 1, c = 0.
    CHECK(brackets(rt::solve<rt::Bisection>(minus_1, rt::Bracket{ -dmax, dmax }), 1.0, 0.0));
    CHECK(brackets(rt::solve<rt::Brent>(minus_1, rt::Bracket{ -dmax, dmax }), 1.0, 0.0));
    CHECK(brackets(rt::solve<rt::Brent>(minus_1, rt::Bracket{ 0.0, dmax }), 1.0, 0.0));

    // A wide bracket of opposite signs whose fn is not finite at 0, so zero first is skipped: the guard's
    // overflow-free midpoint and distances. The root is far from 0, since a method may reach 0 by itself (from
    // [-max, max / 2], Brent's secant towards a root at 1 rounds onto 0; on [-max, max], Bisection's first midpoint
    // is 0), and a NaN there is non_finite like any other. fn changes sign at fl(1e300) exactly: c = 0. The ends
    // stay below max / 2 in magnitude, where fn itself would overflow.
    auto const nan_at_0 = [](double x) { return x == 0.0 ? dnan : x - 1e300; };
    CHECK(brackets(rt::solve<rt::Bisection>(nan_at_0, rt::Bracket{ -dmax / 2.0, dmax / 4.0 }), 1e300, 0.0));
    CHECK(brackets(rt::solve<rt::Brent>(nan_at_0, rt::Bracket{ -dmax / 2.0, dmax / 4.0 }), 1e300, 0.0));
    auto const nan_at_0_root_1 = [](double x) { return x == 0.0 ? dnan : x - 1.0; };
    CHECK(rt::solve<rt::Bisection>(nan_at_0_root_1, rt::Bracket{ -dmax, dmax })
          == std::unexpected(rt::error::non_finite)); // the first midpoint is 0

    // A large root: x - 1e300 is exact near 1e300 in the sense that fn changes sign at fl(1e300).
    CHECK(brackets(rt::solve<rt::Brent>([](double x) { return x - 1e300; }, rt::Bracket{ 0.0, dmax }), 1e300, 0.0));

    // A jump at 0, the documented walk: with fn(0) = 1 the bracket becomes [-1, 0], approached down to denorm_min
    // (1076 evaluations in double, for both methods, since Brent bisects at a jump).
    auto const jump = [](double x) { return x < 0.0 ? -1.0 : 1.0; };
    auto const [jb, njb] = counted<rt::Bisection>(jump, rt::Bracket{ -1.0, 2.0 });
    auto const [jr, njr] = counted<rt::Brent>(jump, rt::Bracket{ -1.0, 2.0 });
    CHECK((jb && nxx::math::abs(jb->value) <= 2.0 * dmin && njb == 1076));
    CHECK((jr && nxx::math::abs(jr->value) <= 2.0 * dmin && njr == 1076));
}

TEST_CASE("roots: poles and steep functions")
{
    // A pole brackets like a root and is returned as one (documented): 1/(x - 1/3) on [0, 1] converges to fl(1/3),
    // where |fn| is huge.
    double const third = 1.0 / 3.0;
    auto const pole = [third](double x) { return 1.0 / (x - third); };
    auto const p = rt::solve<rt::Brent>(pole, rt::Bracket{ 0.0, 1.0 });
    CHECK((p && nxx::math::abs(p->value - third) <= p->abs_error && nxx::math::abs(pole(p->value)) > 1e15));

    // The two continuous functions the rejected pole rule misread: a resonance, whose values near the root exceed its
    // end values (its sign change is exactly that of x^2 - 2, c = 8e-17 as above), and a Gaussian-damped root, whose
    // end values are 2e-49 (x - 1/3 is exact near 1/3, so fn changes sign at fl(1/3), c = 0).
    auto const resonance = [](double x) {
        double const u = x * x - 2.0;
        return u / (u * u + 1e-24);
    };
    CHECK(brackets(rt::solve<rt::Bisection>(resonance, rt::Bracket{ 1.0, 2.0 }), sqrt2, 8e-17));
    CHECK(brackets(rt::solve<rt::Brent>(resonance, rt::Bracket{ 1.0, 2.0 }), sqrt2, 8e-17));
    auto const damped = [third](double x) { return (x - third) * std::exp(-1000.0 * (x - third) * (x - third)); };
    CHECK(brackets(rt::solve<rt::Bisection>(damped, rt::Bracket{ 0.0, 1.0 }), third, 0.0));
    CHECK(brackets(rt::solve<rt::Brent>(damped, rt::Bracket{ 0.0, 1.0 }), third, 0.0));

    // A step function: its switching point, which is the sign change asked for.
    auto const step = [third](double x) { return x < third ? -1.0 : 1.0; };
    auto const s = rt::solve<rt::Brent>(step, rt::Bracket{ 0.0, 1.0 });
    CHECK((s && nxx::math::abs(s->value - third) <= s->abs_error));
}

TEST_CASE("roots: float and long double")
{
    // x^2 - 2 and the Wallis cubic, with fn's rounding (under eps |ref|, as for double above) as c.
    auto const sq_f = [](float x) { return x * x - 2.0f; };
    float const ef = std::numeric_limits<float>::epsilon();
    float const tf = std::numeric_limits<float>::denorm_min();
    auto const rf = rt::solve<rt::Brent>(sq_f, rt::Bracket{ 1.0f, 2.0f });
    CHECK(brackets(rf, 1.41421356f, ef * 1.5f, 2.0f * (tf + 2.0f * ef * 1.5f)));
    auto const bf = rt::solve<rt::Bisection>(sq_f, rt::Bracket{ 1.0f, 2.0f });
    CHECK(brackets(bf, 1.41421356f, ef * 1.5f, 2.0f * (tf + 2.0f * ef * 1.5f)));

    using ld = long double;
    ld const el = std::numeric_limits<ld>::epsilon();
    ld const tl = std::numeric_limits<ld>::denorm_min();
    auto const wl = rt::solve<rt::Brent>([](ld x) { return (x * x - 2.0L) * x - 5.0L; }, rt::Bracket{ 2.0L, 3.0L });
    CHECK(brackets(wl, 2.09455148154232659148238654057930296L, el * 2.1L, 2.0L * (tl + 2.0L * el * 2.1L)));
}

TEST_CASE("roots: multiprecision types")
{
    // x^2 - 2 and the Wallis cubic at 50 digits, against 70-digit references, with c = eps |ref| (fn's rounding, as
    // derived for double: 0.35 eps and 0.41 eps relative).
    auto const check = [](auto zero) {
        using T = decltype(zero);
        T const eps = std::numeric_limits<T>::epsilon();
        T const tiny = T(std::numeric_limits<T>::denorm_min());
        T const s2 = T(sqrt2_70);
        T const w = T(wallis_70);
        auto const sq = [](T const &x) { return T(x * x - 2); };
        auto const cubic = [](T const &x) { return T(T(x * x - 2) * x - 5); };
        for (auto const &r : { rt::solve<rt::Brent>(sq, rt::Bracket<T>{ T(1), T(2) }),
                 rt::solve<rt::Bisection>(sq, rt::Bracket<T>{ T(1), T(2) }) })
            CHECK(brackets(r, s2, T(eps * s2), T(2 * T(tiny + T(2 * eps * s2)))));
        for (auto const &r : { rt::solve<rt::Brent>(cubic, rt::Bracket<T>{ T(2), T(3) }),
                 rt::solve<rt::Bisection>(cubic, rt::Bracket<T>{ T(2), T(3) }) })
            CHECK(brackets(r, w, T(eps * w), T(2 * T(tiny + T(2 * eps * w)))));
        // Criteria from doubles, through type_identity: {.rel_tol = 1e-30}.
        auto const loose = rt::solve<rt::Brent>(sq, rt::Bracket<T>{ T(1), T(2) }, { .rel_tol = 1e-30 });
        CHECK(brackets(loose, s2, T(eps * s2), T(2 * T(tiny + T(T(1e-30) * s2)))));
    };
    check(bin50{});
    check(bin50_et{});
    check(dec50{});

    // Brent's counts for the Wallis cubic at 50 digits are not asserted: no independent model rounds like these types.

    // cpp_bin_float has no subnormals: a difference below its smallest number t flushes to 0, which must not become an
    // abs_error of 0 (the independent review's case, 02fffe5). fn = x / t - q changes sign at t q, between ends
    // t and t (1 + 4 eps) whose difference flushes; x / t is exact (t is a power of 2), so the distance is checked as
    // |value / t - q| <= abs_error / t.
    bin50 const t = std::numeric_limits<bin50>::denorm_min();
    bin50 const e = std::numeric_limits<bin50>::epsilon();
    bin50 const q = bin50(bin50(1) + bin50(2) * e);
    auto const scaled = [&t, &q](bin50 const &x) { return bin50(bin50(x / t) - q); };
    rt::Bracket<bin50> const smallest{ t, bin50(t * bin50(bin50(1) + bin50(4) * e)) };
    for (auto const &r : { rt::solve<rt::Bisection>(scaled, smallest), rt::solve<rt::Brent>(scaled, smallest) }) {
        REQUIRE(r.has_value());
        CHECK(r->abs_error > bin50(0));
        CHECK(bin50(abs(bin50(bin50(r->value / t) - q))) <= bin50(r->abs_error / t));
    }
}

// Emscripten's default build compiles throw but aborts on it rather than unwinding to a catch, so the test runs
// everywhere else.
#ifndef __EMSCRIPTEN__
TEST_CASE("roots: fn's exceptions pass through")
{
    auto const throwing = [](double x) -> double {
        if (x > 1.2) throw std::runtime_error("fn failed");
        return x * x - 2.0;
    };
    CHECK_THROWS_AS((void)rt::solve<rt::Brent>(throwing, rt::Bracket{ 1.0, 2.0 }), std::runtime_error);
}
#endif

TEST_CASE("roots: error formats as its name")
{
    CHECK(std::format("{}", rt::error::non_finite) == "non_finite");
    CHECK(std::format("{}", rt::error::no_sign_change) == "no_sign_change");
    CHECK(std::format("[{:>16}]", rt::error::no_sign_change) == "[  no_sign_change]");
}
