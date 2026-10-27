// Ridders, RegulaFalsi, Illinois and Steffensen in nxx::roots (Phase 6, docs/roots/PHASE6_API.md), and the square root
// Ridders needs.
//
// Evaluation counts are asserted where no libm is involved (+ - * / and sqrt in double), against an independent Python
// implementation of the specification (.scratch/roots6/proto6.py, with detail::sqrt ported step for step), assuming
// no FMA contraction, as in Phases 2 to 5. Every other case checks accuracy (roots_check.hpp).

#include "roots_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <expected>
#include <limits>
#include <tuple>

using namespace roots_check;

namespace {

    constexpr auto sq = [](double x) { return x * x - 2.0; };
    constexpr auto wal = [](double x) { return (x * x - 2.0) * x - 5.0; };
    constexpr auto x20 = [](double x) {
        double const x2 = x * x;
        double const x4 = x2 * x2;
        double const x8 = x4 * x4;
        return x8 * x8 * x4 - 1.0;
    };
    constexpr auto triple = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };
    constexpr auto tiny_root = [](double x) { return x - 1e-300; };
    constexpr auto minus_1 = [](double x) { return x - 1.0; };
    constexpr auto cycle = [](double x) { return (x * x - 2.0) * x + 2.0; };
    constexpr auto reciprocal = [](double x) { return 1.0 / x - 1.0; };
    constexpr auto jump = [](double x) { return x < 0.0 ? -1.0 : 1.0; };
    constexpr auto pole = [](double x) { return 1.0 / (x - 1.0 / 3.0); };
    auto const sqrt_fn = [](double x) { return std::sqrt(x) - 0.5; };
    constexpr auto w10 = [](double x) {
        constexpr double c[] = { 1.0, -55.0, 1320.0, -18150.0, 157773.0, -902055.0, 3416930.0, -8409500.0, 12753576.0,
            -10628640.0, 3628800.0 };
        double p = 0.0;
        for (double const a : c) p = p * x + a;
        return p;
    };

    constexpr double cycle_root = -1.769292354238631415240409; // x^3 - 2x + 2 (mpmath)

    // y^n by repeated multiplication, as the prototype; Alefeld, Potra and Shi's problems (Table I).
    constexpr double power(double y, int n)
    {
        double r = y;
        for (int i = 1; i < n; ++i) r = r * y;
        return r;
    }
    constexpr double problem2(double x)
    {
        double s = 0.0;
        for (int i = 1; i <= 20; ++i) {
            double const t = x - double(i * i);
            s = s + double((2 * i - 5) * (2 * i - 5)) / (t * t * t);
        }
        return -2.0 * s;
    }
    auto problem4(int n) { return [n](double x) { return power(x, n) - 1.0; }; }
    auto problem7(int n)
    {
        return [n](double x) {
            double const m = double(n);
            return double(1 + (1 - n) * (1 - n)) * x - (1.0 - m * x) * (1.0 - m * x);
        };
    }
    auto problem8(int n) { return [n](double x) { return x * x - power(1.0 - x, n); }; }
    auto problem9(int n)
    {
        return [n](double x) { return (1.0 + power(double(1 - n), 4)) * x - power(1.0 - double(n) * x, 4); };
    }
    auto problem11(int n) { return [n](double x) { return (double(n) * x - 1.0) / (double(n - 1) * x); }; }

    template<typename M, typename S>
    concept solvable = requires(double (*f)(double), S s) { rt::solve<M>(f, s); };

    using B = rt::Bracket<double>;

} // namespace

// The square root (detail::sqrt): exact on perfect squares and even powers of 2, and 0, +inf as they are.
static_assert(rt::detail::sqrt(4.0) == 2.0 && rt::detail::sqrt(144.0) == 12.0 && rt::detail::sqrt(0.25) == 0.5);
static_assert(rt::detail::sqrt(0x1p-1074) == 0x1p-537 && rt::detail::sqrt(0x1p1022) == 0x1p511);
static_assert(rt::detail::sqrt(0.0) == 0.0 && rt::detail::sqrt(dinf) == dinf);
static_assert(rt::detail::sqrt(9.0f) == 3.0f && rt::detail::sqrt(2.25L) == 1.5L);

// Constant evaluation: Ridders and Steffensen in each type, RegulaFalsi and Illinois in double (their loop's arithmetic
// is AndersonBjorck's, tested in each type). The bound is the contract's, 2 tol + c, about 7 eps near sqrt 2.
static_assert(nxx::math::abs(rt::solve<rt::Ridders>(sq, B{ 1.0, 2.0 })->value - sqrt2) <= 8 * deps);
static_assert(nxx::math::abs(rt::solve<rt::Steffensen>(sq, 1.0)->value - sqrt2) <= 8 * deps);
static_assert(nxx::math::abs(rt::solve<rt::RegulaFalsi>(sq, B{ 1.0, 2.0 })->value - sqrt2) <= 8 * deps);
static_assert(nxx::math::abs(rt::solve<rt::Illinois>(sq, B{ 1.0, 2.0 })->value - sqrt2) <= 8 * deps);
static_assert(nxx::math::abs(
                  rt::solve<rt::Ridders>([](float x) { return x * x - 2.0f; }, rt::Bracket{ 1.0f, 2.0f })->value
                  - 1.41421356f)
              <= 8 * std::numeric_limits<float>::epsilon());
static_assert(nxx::math::abs(rt::solve<rt::Steffensen>([](float x) { return x * x - 2.0f; }, 1.0f)->value - 1.41421356f)
              <= 8 * std::numeric_limits<float>::epsilon());
static_assert(nxx::math::abs(rt::solve<rt::Ridders>([](long double x) { return x * x - 2.0L; },
                                 rt::Bracket{ 1.0L, 2.0L })
                                 ->value
                             - 1.4142135623730950488016887242096980785696718753769L)
              <= 8 * std::numeric_limits<long double>::epsilon());
static_assert(nxx::math::abs(rt::solve<rt::Steffensen>([](long double x) { return x * x - 2.0L; }, 1.0L)->value
                             - 1.4142135623730950488016887242096980785696718753769L)
              <= 8 * std::numeric_limits<long double>::epsilon());

// Each new tag at the wrong start is rejected.
static_assert(solvable<rt::Ridders, B> && !solvable<rt::Ridders, double>);
static_assert(solvable<rt::RegulaFalsi, B> && !solvable<rt::RegulaFalsi, double>);
static_assert(solvable<rt::Illinois, B> && !solvable<rt::Illinois, double>);
static_assert(solvable<rt::Steffensen, double> && !solvable<rt::Steffensen, B>);

TEST_CASE("roots: Ridders, Illinois and RegulaFalsi, evaluation counts in double")
{
    auto const check = [](auto const &fn, B const &bracket, int ridders, int illinois, int regula_falsi) {
        CHECK(counted<rt::Ridders>(fn, bracket).second == ridders);
        CHECK(counted<rt::Illinois>(fn, bracket).second == illinois);
        CHECK(counted<rt::RegulaFalsi>(fn, bracket).second == regula_falsi);
    };
    check(sq, { 1.0, 2.0 }, 12, 11, 20);
    check(wal, { 2.0, 3.0 }, 14, 10, 23);
    check(x20, { 0.0, 5.0 }, 18, 30, 45); // RegulaFalsi's guard: over 200,000 evaluations without it
    check(triple, { 0.0, 3.0 }, 86, 148, 145); // Ridders on a triple root; the guard makes RegulaFalsi cheaper here
    check(tiny_root, { 0.0, 1.0 }, 4, 4, 4); // Ridders' point without cancellation: 1922 in the textbook form
    check(minus_1, { -dmax, dmax }, 9, 13, 13);
    check(minus_1, { 0.0, dmax }, 8, 12, 12);
    check(cycle, { -2.0, -1.0 }, 14, 12, 20);
    check(reciprocal, { 0.5, 3.0 }, 12, 12, 25);
    check(w10, { 6.7, 7.6 }, 18, 24, 14);
    check(jump, { -1.0, 2.0 }, 778, 51, 1076); // the walk to 0
    check(pole, { 0.0, 1.0 }, 79, 95, 197);
    check(sqrt_fn, { 0.0, 1.0 }, 12, 12, 22);

    // Alefeld, Potra and Shi's problems.
    for (auto const &[n, r, i, f] : { std::tuple{ 1, 18, 19, 26 }, std::tuple{ 5, 18, 14, 17 }, std::tuple{ 10, 16, 17, 21 } })
        check(problem2, { n * n + 1e-9, (n + 1) * (n + 1) - 1e-9 }, r, i, f);
    for (auto const &[n, r, i, f] : { std::tuple{ 4, 16, 21, 34 }, std::tuple{ 8, 14, 23, 37 }, std::tuple{ 12, 16, 24, 37 } })
        check(problem4(n), { 0.0, 5.0 }, r, i, f);
    for (auto const &[n, r, i, f] : { std::tuple{ 5, 18, 10, 17 }, std::tuple{ 10, 18, 9, 14 }, std::tuple{ 20, 18, 9, 12 } })
        check(problem7(n), { 0.0, 1.0 }, r, i, f);
    for (auto const &[n, r, i, f] : { std::tuple{ 5, 14, 10, 24 }, std::tuple{ 10, 16, 13, 26 }, std::tuple{ 20, 14, 12, 22 } })
        check(problem8(n), { 0.0, 1.0 }, r, i, f);
    for (auto const &[n, r, i, f] : { std::tuple{ 5, 12, 9, 11 }, std::tuple{ 15, 14, 9, 8 }, std::tuple{ 20, 14, 9, 8 } })
        check(problem9(n), { 0.0, 1.0 }, r, i, f);
    for (auto const &[n, r, i, f] : { std::tuple{ 5, 16, 15, 35 }, std::tuple{ 15, 16, 17, 32 }, std::tuple{ 20, 18, 16, 33 } })
        check(problem11(n), { 0.01, 1.0 }, r, i, f);
}

TEST_CASE("roots: Ridders, Illinois and RegulaFalsi, accuracy in double")
{
    auto const all = [](auto const &fn, B const &bracket, double ref, double c) {
        CHECK(brackets(rt::solve<rt::Ridders>(fn, bracket), ref, c));
        CHECK(brackets(rt::solve<rt::Illinois>(fn, bracket), ref, c));
        CHECK(brackets(rt::solve<rt::RegulaFalsi>(fn, bracket), ref, c));
    };
    // The conditioning terms of Phases 1 to 4: x^2 - 2 7.9e-17; Wallis 2e-16; x^20 - 1 1.1e-16; x - 1, (x - 1)^3 and
    // x - 1e-300 exact; x^3 - 2x + 2 1.8e-16; 1/x - 1 1.2e-16; sqrt exactly 0 at 0.25.
    all(sq, { 1.0, 2.0 }, sqrt2, 8e-17);
    all(wal, { 2.0, 3.0 }, wallis, 2e-16);
    all(x20, { 0.0, 5.0 }, 1.0, 1.1e-16);
    all(triple, { 0.0, 3.0 }, 1.0, 0.0);
    all(tiny_root, { 0.0, 1.0 }, 1e-300, 0.0);
    all(minus_1, { -dmax, dmax }, 1.0, 0.0);
    all(cycle, { -2.0, -1.0 }, cycle_root, 1.8e-16);
    all(reciprocal, { 0.5, 3.0 }, 1.0, 1.2e-16);
    all(sqrt_fn, { 0.0, 1.0 }, 0.25, 0.0);
    // libm, with Phase 1's budget of 2 ulps.
    all([](double x) { return std::cos(x) - x; }, { 0.0, 1.0 }, dottie, 1.4e-16);
    all([](double x) { return std::exp(x) - 2.0; }, { 0.0, 1.0 }, ln2, 4.5e-16);
    all([](double x) { return x * std::exp(x) - 1.0; }, { 0.0, 1.0 }, omega, 1.3e-16);
    all([](double x) { return std::log(x / 1e-12); }, { 1e-20, 1.0 }, 1e-12, 2.5e-28);
}

TEST_CASE("roots: Steffensen")
{
    // The prototype's counts and the root found, with Secant's for contrast: f's scale changes the cost and the root.
    auto const check = [](auto const &fn, double x0, int count, double root) {
        auto const [r, n] = counted<rt::Steffensen>(fn, x0);
        REQUIRE(r.has_value());
        CHECK(n == count);
        CHECK(nxx::math::abs(r->value - root) <= r->abs_error + 2e-16); // the conditioning terms above
    };
    check(sq, 1.0, 10, sqrt2);
    check(sq, 10.0, 39, sqrt2); // from far: Secant 13
    check(sq, 0.5, 11, -sqrt2); // the first step lands at -1.83; Secant and Newton find sqrt 2
    check(wal, 2.0, 9, wallis);
    check(wal, 3.0, 45, wallis);
    check([](double x) { return 1000.0 * (x * x - 2.0); }, 1.0, 27, -sqrt2); // z = 1 - 1000 brackets -sqrt 2
    check([](double x) { return (x * x - 2.0) / 1000.0; }, 1.0, 9, sqrt2);
    check(reciprocal, 0.5, 10, 1.0);
    auto const [past_pole, n] = counted<rt::Steffensen>(reciprocal, 3.0);
    CHECK(past_pole == std::unexpected(rt::error::not_converged));
    CHECK(n == 12);
    CHECK(counted<rt::Secant>([](double x) { return 1000.0 * (x * x - 2.0); }, 1.0).first->value > 0.0);

    // f(x) below the resolution of x: z keeps a distance tol(x), and the root is found.
    CHECK(brackets(rt::solve<rt::Steffensen>(sq, sqrt2), sqrt2, 8e-17));
    // x + f(x) would overflow: not_converged after one evaluation.
    auto const [overflow, n_overflow] = counted<rt::Steffensen>([](double) { return 1e308; }, 1e308);
    CHECK(overflow == std::unexpected(rt::error::not_converged));
    CHECK(n_overflow == 1);
    // fn not finite at z (z leaves fn's domain): non_finite.
    CHECK(rt::solve<rt::Steffensen>([](double x) { return std::log(x) - 5.0; }, 0.5)
          == std::unexpected(rt::error::non_finite));
    // libm: cos x - x from 0 (9 evaluations in the prototype; accuracy only).
    CHECK(brackets(rt::solve<rt::Steffensen>([](double x) { return std::cos(x) - x; }, 0.0), dottie, 1.4e-16));
}

TEST_CASE("roots: the square root and the new methods in other types")
{
    // detail::sqrt within one ulp of std::sqrt (correctly rounded by IEEE 754), over a range of magnitudes.
    for (double const x : { 2.0, 3.0, 0.1, 1e-300, 1e300, dmin, dmax, 0x1.fffffffffffffp-1, 7.0 }) {
        double const a = rt::detail::sqrt(x);
        double const b = std::sqrt(x);
        CHECK(std::abs(a - b) <= std::nextafter(b, dinf) - b);
    }
    for (float const x : { 2.0f, 3.0f, 0.1f, 1e-38f, 1e38f, std::numeric_limits<float>::denorm_min() }) {
        float const a = rt::detail::sqrt(x);
        float const b = std::sqrt(x);
        CHECK(std::abs(a - b) <= std::nextafter(b, std::numeric_limits<float>::infinity()) - b);
    }
    // cpp_bin_float_50 at 2 and at 2e-1000000 (the logarithmic scaling), cpp_dec_float_50 at 2: within one ulp of the
    // 70-digit sqrt 2 (scaled).
    bin50 const s2 = from_string<bin50>(sqrt2_70);
    bin50 const bin_eps = std::numeric_limits<bin50>::epsilon();
    CHECK(abs(bin50(rt::detail::sqrt(bin50(2)) - s2)) <= bin50(bin_eps * s2));
    bin50 const tiny = rt::detail::sqrt(bin50("2e-1000000"));
    CHECK(abs(bin50(tiny / bin50("1e-500000") - s2)) <= bin50(2 * bin_eps * s2));
    dec50 const d2 = from_string<dec50>(sqrt2_70);
    CHECK(abs(dec50(rt::detail::sqrt(dec50(2)) - d2)) <= dec50(std::numeric_limits<dec50>::epsilon() * d2));

    // The methods: x^2 - 2 with c = eps |ref| (fn's rounding), as for Brent.
    auto const check = [](auto zero) {
        using T = decltype(zero);
        T const eps = std::numeric_limits<T>::epsilon();
        T const t = T(std::numeric_limits<T>::denorm_min());
        T const ref = from_string<T>(sqrt2_70);
        auto const sq_t = [](T const &x) { return T(x * x - 2); };
        T const bound = T(2 * T(t + T(2 * eps * ref)));
        CHECK(brackets(rt::solve<rt::Ridders>(sq_t, rt::Bracket<T>{ T(1), T(2) }), ref, T(eps * ref), bound));
        CHECK(brackets(rt::solve<rt::Steffensen>(sq_t, T(1)), ref, T(eps * ref), bound));
    };
    check(bin50{});
    check(bin50_et{});
    check(dec50{});
    auto const sq_et = [](bin50_et const &x) { return bin50_et(x * x - 2); };
    bin50_et const ref_et = from_string<bin50_et>(sqrt2_70);
    bin50_et const eps_et = std::numeric_limits<bin50_et>::epsilon();
    for (auto const &r : { rt::solve<rt::RegulaFalsi>(sq_et, rt::Bracket<bin50_et>{ 1, 2 }),
             rt::solve<rt::Illinois>(sq_et, rt::Bracket<bin50_et>{ 1, 2 }) })
        CHECK(brackets(r, ref_et, bin50_et(eps_et * ref_et),
            bin50_et(2 * bin50_et(bin50_et(std::numeric_limits<bin50_et>::denorm_min()) + bin50_et(2 * eps_et * ref_et)))));
}
