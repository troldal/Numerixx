// Ridders' extrapolation (derivative<Ridders>) and the Estimate type every method returns.
//
// Accuracy bounds (docs/PHASE3_API.md, "Tests"): the actual error must be within 6 times the estimate (the Python
// prototype's worst over 5,400 random cases was 3.6) and within 4 times the prototype's actual error for the same case
// (4 eps where that was 0); the estimate within twice the prototype's; the evaluations at most the prototype's plus 4
// (libm implementations round fn's values differently, which can move a stall by a level). The prototype is an
// independent implementation; its values are in the table of docs/PHASE3_API.md. For the multiprecision types, which
// the prototype does not model, the estimate is bounded by twice its value in the first C++ run.

#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <doctest/doctest.h>

#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <type_traits>

namespace mp = boost::multiprecision;

using bin50 = mp::cpp_bin_float_50;
using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
using dec50 = mp::cpp_dec_float_50;

using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Estimate;
using nxx::deriv::Ridders;

namespace {

    constexpr auto square = [](auto x) { return x * x; };
    constexpr auto cube = [](auto x) { return x * x * x; };
    constexpr auto identity = [](double x) { return x; };
    constexpr auto huge = [](double) { return 1e300; };
    constexpr auto not_a_number = [](double) { return std::numeric_limits<double>::quiet_NaN(); };
    constexpr double dmax = std::numeric_limits<double>::max();
    constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr auto non_finite = std::unexpected(nxx::deriv::error::non_finite);

    // Ridders' rounding floor at step h, by its definition: (eps |f(x - h)| + eps |f(x + h)|) / ((x + h) - (x - h)).
    template<typename T, typename F>
    constexpr T floor_at(F f, T x, T h)
    {
        T const eps = std::numeric_limits<T>::epsilon();
        T const lo = x - h;
        T const hi = x + h;
        return (eps * nxx::math::abs(T(f(lo))) + eps * nxx::math::abs(T(f(hi)))) / (hi - lo);
    }

    // The prototype's values for a case, relative to |f'|, and its number of evaluations.
    struct Reference
    {
        double error;
        double estimate;
        int evaluations;
    };

    // fn wrapped to count its calls.
    template<typename F>
    struct Counted
    {
        F fn;
        int *calls;
        template<typename T>
        constexpr auto operator()(T const &x) const
        {
            ++*calls;
            return fn(x);
        }
    };

    // The accuracy check of the header comment, for double and float.
    template<typename T, typename F>
    bool accurate(F fn, T x, T exact, Reference ref)
    {
        int calls = 0;
        auto const r = derivative<Ridders>(Counted{ fn, &calls }, x);
        if (!r || !r->abs_error) return false;
        T const eps = std::numeric_limits<T>::epsilon();
        T const actual = std::abs(r->value - exact);
        T const error_bound = ref.error > 0 ? T(4) * T(ref.error) * std::abs(exact) : T(4) * eps * std::abs(exact);
        bool const ok = actual <= T(6) * *r->abs_error && actual <= error_bound
                        && *r->abs_error <= T(2) * T(ref.estimate) * std::abs(exact) && calls <= ref.evaluations + 4;
        if (!ok)
            MESSAGE("actual ", actual, ", estimate ", *r->abs_error, ", error bound ", error_bound, ", calls ", calls);
        return ok;
    }

    double some_function(double x) { return x + 1; }

    // f = (max h) (a t + b t^3) with t = (x - 1) / h, so f'(1) = a max, and the central difference at step s is
    // D(s) = max (a + b s^2 / h^2): every tableau entry is finite, but two of opposite sign can differ by more than
    // max. a = 0.9, b = -1.2: D(h) = -0.3 max and D(h/2) = 0.6 max extrapolate to A[1][1] = 0.9 max, which differs
    // from its parent D(h) by 1.2 max.
    constexpr double h_change = 0x1p-51;
    constexpr auto opposite_change = [](double x) {
        double const t = (x - 1.0) / h_change;
        return (dmax * h_change) * (0.9 * t - 1.2 * t * t * t);
    };
    // a = 0.5, b = -1.4: D(h) = -0.9 max and D(h/2) = 0.15 max differ by 1.05 max within the tableau update.
    constexpr double h_levels = 0x1p-20;
    constexpr auto opposite_levels = [](double x) {
        double const t = (x - 1.0) / h_levels;
        return (dmax * h_levels) * (0.5 * t - 1.4 * t * t * t);
    };
    constexpr double deps = std::numeric_limits<double>::epsilon();

    constexpr bool within(double a, double b, double tolerance) { return nxx::math::abs(a - b) <= tolerance; }

} // namespace

// Polynomials with dyadic x and h are exact at every level (binary types only: decimal division need not round
// correctly). For x^3 at 2 with h = 1/16, A[2][2] equals its parent A[1][1] = 12 exactly, so the tableau stops there
// with the rounding floor of level 2 (h = 1/64) as the estimate. For x^2, A[1][1] = A[0][0] = 4 at level 1 (h = 1/32).
static_assert(derivative<Ridders>(cube, 2.0, 1.0 / 16)->value == 12.0);
static_assert(*derivative<Ridders>(cube, 2.0, 1.0 / 16)->abs_error == floor_at(cube, 2.0, 1.0 / 64));
static_assert(derivative<Ridders>(square, 2.0, 1.0 / 16)->value == 4.0);
static_assert(*derivative<Ridders>(square, 2.0, 1.0 / 16)->abs_error == floor_at(square, 2.0, 1.0 / 32));
// float: x = 1/2 with h = 1/4 keeps every point and value exact in 24 bits.
static_assert(derivative<Ridders>(cube, 0.5f, 0.25f)->value == 0.75f);
static_assert(*derivative<Ridders>(cube, 0.5f, 0.25f)->abs_error == floor_at(cube, 0.5f, 1.0f / 16));

// Failures, also in constant expressions.
static_assert(derivative<Ridders>(identity, qnan) == non_finite);
static_assert(derivative<Ridders>(identity, inf) == non_finite);
static_assert(derivative<Ridders>(identity, -inf, 1.0) == non_finite);
static_assert(derivative<Ridders>(not_a_number, 1.0) == non_finite); // never finite: 32 failing levels
static_assert(derivative<Ridders>(identity, dmax) == non_finite); // every step fails the overflow check
static_assert(derivative<Ridders>(identity, 1.0, 1e-300) == non_finite); // h below the resolution of x at once

// Regression: the rounding floor of a large fn over a tiny step, eps 2e300 / 2e-300, overflows; it saturates at max,
// and is computed for Ridders only, so the stencil result is a constant expression again.
static_assert(derivative(huge, 0.0, 1e-300) == 0.0);
static_assert(derivative<Ridders>(huge, 0.0, 1e-300) == Estimate<double>{ .value = 0.0, .abs_error = dmax });

// Regression: finite tableau entries of opposite sign near max differ by more than max. The change from the parent
// and the disagreement overflowed to an infinite estimate (a compile error here), and the tableau update overflowed
// (a compile error here, a restart at run time). Now the change saturates at max, like the floor: the result is
// untrusted, the restart below h_change meets the resolution of x, and the first tableau's result is returned. Its
// value is 0.9 max to within the rounding of fn's few operations, a few eps relative, so 4 eps.
static_assert(derivative<Ridders>(opposite_change, 1.0, h_change)->abs_error == dmax);
static_assert(within(derivative<Ridders>(opposite_change, 1.0, h_change)->value, 0.9 * dmax, 4 * deps * dmax));
// The tableau update that would exceed max counts as a non-finite level, so the tableau restarts below it at h/4,
// where the cubic is exact up to rounding: the estimate covers the actual error and is near the rounding floor, about
// eps relative to max / 2.
static_assert(within(derivative<Ridders>(opposite_levels, 1.0, h_levels)->value,
    dmax / 2,
    *derivative<Ridders>(opposite_levels, 1.0, h_levels)->abs_error));
static_assert(*derivative<Ridders>(opposite_levels, 1.0, h_levels)->abs_error <= 4 * deps * (dmax / 2));

// The uniform return type: stencils leave the error estimate empty, Ridders fills it in.
static_assert(!derivative(cube, 2.0, 0.25)->abs_error.has_value());
static_assert(!derivative<Central<4>>(cube, 2.0, 0.25)->abs_error.has_value());
static_assert(derivative<Ridders>(cube, 2.0, 1.0 / 16)->abs_error.has_value());

// Estimate: an aggregate, convertible to T, compared member by member. The cases of the conversion spike.
static_assert(std::is_aggregate_v<Estimate<double>>);
namespace {
    constexpr Estimate<double> one{ .value = 1.0, .abs_error = std::nullopt };
    constexpr double converted = *derivative(cube, 2.0, 0.25);
    constexpr double twice(double x) { return 2 * x; }
} // namespace
static_assert(one == 1.0 && one + 1.0 == 2.0 && converted == 12.0625 && twice(one) == 2.0);
static_assert(derivative(cube, 2.0, 0.25) == 12.0625);
static_assert(*derivative<Ridders>(cube, 2.0, 1.0 / 16) == 12.0);
static_assert(derivative<Ridders>(cube, 2.0, 1.0 / 16).transform([](double d) { return d / 2; }) == 6.0);
static_assert(one != Estimate<double>{ .value = 1.0, .abs_error = 0.5 }); // the error is compared too

TEST_CASE("ridders: polynomials, the uniform type and the conversion at run time")
{
    // The compile-time cases above.
    CHECK(derivative<Ridders>(cube, 2.0, 1.0 / 16)->value == 12.0);
    CHECK(*derivative<Ridders>(cube, 2.0, 1.0 / 16)->abs_error == floor_at(cube, 2.0, 1.0 / 64));
    CHECK(derivative<Ridders>(square, 2.0, 1.0 / 16)->value == 4.0);
    CHECK(*derivative<Ridders>(square, 2.0, 1.0 / 16)->abs_error == floor_at(square, 2.0, 1.0 / 32));
    CHECK(derivative<Ridders>(cube, 0.5f, 0.25f)->value == 0.75f);
    CHECK(*derivative<Ridders>(cube, 0.5f, 0.25f)->abs_error == floor_at(cube, 0.5f, 1.0f / 16));
    CHECK(derivative(huge, 0.0, 1e-300) == 0.0);
    CHECK(derivative<Ridders>(huge, 0.0, 1e-300) == Estimate<double>{ .value = 0.0, .abs_error = dmax });
    auto const change = derivative<Ridders>(opposite_change, 1.0, h_change);
    REQUIRE(change.has_value());
    CHECK(change->abs_error == dmax);
    CHECK(within(change->value, 0.9 * dmax, 4 * deps * dmax));
    auto const levels = derivative<Ridders>(opposite_levels, 1.0, h_levels);
    REQUIRE(levels.has_value());
    CHECK(within(levels->value, dmax / 2, *levels->abs_error));
    CHECK(*levels->abs_error <= 4 * deps * (dmax / 2));
    // Compile time and run time agree.
    constexpr auto levels_at_compile_time = derivative<Ridders>(opposite_levels, 1.0, h_levels);
    CHECK(levels == levels_at_compile_time);
    CHECK(!derivative(cube, 2.0, 0.25)->abs_error.has_value());
    CHECK(!derivative<Central<4>>(cube, 2.0, 0.25)->abs_error.has_value());
    CHECK(derivative<Ridders>(cube, 2.0, 1.0 / 16)->abs_error.has_value());

    // The exact stop costs 3 levels for x^3 and 2 for x^2.
    int calls = 0;
    CHECK(derivative<Ridders>(Counted{ cube, &calls }, 2.0, 1.0 / 16).has_value());
    CHECK(calls == 6);
    calls = 0;
    CHECK(derivative<Ridders>(Counted{ square, &calls }, 2.0, 1.0 / 16).has_value());
    CHECK(calls == 4);

    // The conversion and structured bindings.
    auto const [value, abs_error] = *derivative<Ridders>(cube, 2.0, 1.0 / 16);
    CHECK(value == 12.0);
    CHECK(abs_error.has_value());
    double const d = *derivative(cube, 2.0, 0.25);
    CHECK(d == 12.0625);
    CHECK(one + 1.0 == 2.0);
    CHECK(twice(one) == 2.0);
    CHECK(some_function(*derivative<Ridders>(cube, 2.0, 1.0 / 16)) == 13.0);
    CHECK(derivative<Ridders>(cube, 2.0, 1.0 / 16).transform([](double v) { return v / 2; }) == 6.0);
    CHECK(one != Estimate<double>{ .value = 1.0, .abs_error = 0.5 });
}

TEST_CASE("ridders: accuracy against closed forms, double and float")
{
    // The prototype's values for the same cases, with the default initial step (docs/PHASE3_API.md).
    auto const check = []<typename T>(T, Reference exp_ref, Reference sin_ref, Reference runge_ref, Reference log_ref) {
        using std::cos, std::exp, std::log, std::sin;
        CHECK(accurate([](T x) { return T(exp(x)); }, T(1), T(exp(T(1))), exp_ref));
        CHECK(accurate([](T x) { return T(sin(x)); }, T(0.5), T(cos(T(0.5))), sin_ref));
        T const x = T(0.3);
        T const runge_exact = T(-50) * x / ((T(1) + T(25) * x * x) * (T(1) + T(25) * x * x));
        CHECK(accurate([](T t) { return T(T(1) / (T(1) + T(25) * t * t)); }, x, runge_exact, runge_ref));
        CHECK(accurate([](T t) { return T(log(t)); }, T(0.2), T(T(1) / T(0.2)), log_ref));
    };
    check(double{},
        { 9.8e-16, 6.0e-14, 12 },
        { 8.7e-15, 3.1e-14, 12 },
        { 5.6e-15, 4.9e-14, 16 },
        { 1.6e-14, 7.3e-14, 16 });
    check(float{}, { 7.2e-6, 1.9e-5, 10 }, { 4.4e-6, 3.3e-5, 10 }, { 3.3e-7, 5.7e-6, 10 }, { 1.1e-6, 9.8e-6, 12 });
}

TEST_CASE("ridders: restarts below a failing or untrusted step")
{
    // Domain edges: the default initial step 1/32 reaches below 0, where log is NaN.
    auto const log = [](double x) { return std::log(x); };
    CHECK(accurate(log, 0.01, 100.0, { 5.0e-14, 1.8e-13, 24 }));
    CHECK(accurate(log, 1e-4, 1e4, { 1.6e-14, 2.8e-13, 36 }));

    // A pole on a sample (x + h = 1/128 at level 2) and one between the samples (non-dyadic, no infinite value).
    CHECK(accurate([](double x) { return 1 / (x - 1.0 / 128); }, 0.0, -16384.0, { 5.3e-15, 5.7e-14, 26 }));
    CHECK(accurate([](double x) { return 1 / (x - 0.01); }, 0.0, -1e4, { 3.8e-15, 3.9e-14, 30 }));

    // Aliasing: sin(64 pi x) vanishes at every sample of the first tableau (h = 1/32, 1/64), which looks like a
    // derivative of 0; the untrusted-result test restarts below it.
    double const omega = 64 * std::numbers::pi;
    CHECK(accurate([omega](double x) { return std::sin(omega * x); }, 0.0, omega, { 0.0, 2.2e-16, 24 }));

    // Finite values whose difference overflows: at x = 1e300 the first steps give D_i = inf, and the tableau restarts.
    // fn is linear, so it has no truncation error and the estimate is the rounding floor, eps |f'| (the prototype was
    // not run on this case).
    int calls = 0;
    auto const steep = derivative<Ridders>(Counted{ [](double x) { return 5e9 * (x - 1e300); }, &calls }, 1e300);
    double const eps = std::numeric_limits<double>::epsilon();
    REQUIRE(steep.has_value());
    CHECK(std::abs(steep->value - 5e9) <= *steep->abs_error);
    CHECK(*steep->abs_error <= 4 * eps * 5e9);
    CHECK(calls <= 64);
}

TEST_CASE("ridders: failures, untrusted results and the level cap")
{
    int calls = 0;
    auto const counted = [&calls](auto f) {
        calls = 0;
        return Counted{ f, &calls };
    };

    // The constant-expression failures above.
    CHECK(derivative<Ridders>(identity, qnan) == non_finite);
    CHECK(derivative<Ridders>(identity, inf) == non_finite);
    CHECK(derivative<Ridders>(identity, -inf, 1.0) == non_finite);

    // fn never finite: every level fails, and the cap ends it after 32 levels.
    CHECK(derivative<Ridders>(counted(not_a_number), 1.0) == non_finite);
    CHECK(calls == 64);
    // x = max: every step fails the overflow check, so fn is never called.
    CHECK(derivative<Ridders>(counted(identity), dmax) == non_finite);
    CHECK(calls == 0);
    // h below the resolution of x, and a non-finite x: no evaluation, no assert.
    CHECK(derivative<Ridders>(counted(identity), 1.0, 1e-300) == non_finite);
    CHECK(calls == 0);
    CHECK(derivative<Ridders>(counted(identity), qnan) == non_finite);
    CHECK(calls == 0);

    // A step at the resolution of x: the Heaviside step at 1 with h = 2^-51 gives D_0 = 2^50, D_1 = 2^51, an untrusted
    // result (estimate above |value| / 4), but a restart cannot go lower, so it is returned after 2 levels.
    auto const step_at_1 = [](double x) { return x < 1.0 ? 0.0 : 1.0; };
    auto const at_resolution = derivative<Ridders>(counted(step_at_1), 1.0, 0x1p-51);
    REQUIRE(at_resolution.has_value());
    CHECK(*at_resolution->abs_error >= std::abs(at_resolution->value) / 4);
    CHECK(calls == 4);

    // The level cap with untrusted results: the Heaviside step at 0 is untrusted at every restart (each D_i doubles),
    // so the cap ends it after 32 levels, and the last result is returned with an estimate above |value| / 4.
    auto const step_at_0 = [](double x) { return x < 0.0 ? 0.0 : 1.0; };
    auto const capped = derivative<Ridders>(counted(step_at_0), 0.0);
    REQUIRE(capped.has_value());
    CHECK(4 * *capped->abs_error >= std::abs(capped->value));
    CHECK(calls == 64);
}

TEST_CASE("ridders: multiprecision types")
{
    // The estimate's ceiling is a regression baseline (the test rule's exception in CLAUDE.md): the relative estimates
    // of the first C++ run, commit 1023248: bin50 exp 4.4e-47, log 3.7e-46; dec50 exp 4.1e-46, log 2.1e-45. The
    // prototype does not model these types; cpp_dec_float's guard digits keep its rounding below its epsilon, so its
    // floor is conservative. The coverage check (actual within 6 times the estimate) uses the closed forms.
    auto const check = []<typename T>(T, char const *exp_reference, char const *log_reference) {
        using std::abs;
        auto const accurate_mp = [](auto const &r, T const &exact, T const &reference) {
            if (!r || !r->abs_error) return false;
            T const actual = abs(T(r->value - exact));
            return actual <= T(6) * *r->abs_error && *r->abs_error <= T(T(2) * reference * abs(exact));
        };
        auto const exp = [](T const &x) { return T(mp::exp(x)); };
        auto const log = [](T const &x) { return T(mp::log(x)); };
        CHECK(accurate_mp(derivative<Ridders>(exp, T(1)), T(mp::exp(T(1))), T(exp_reference)));
        CHECK(accurate_mp(derivative<Ridders>(log, T("0.2")), T(5), T(log_reference)));
    };
    check(bin50{}, "4.4e-47", "3.7e-46");
    check(bin50_et{}, "4.4e-47", "3.7e-46");
    check(dec50{}, "4.1e-46", "2.1e-45");
}
