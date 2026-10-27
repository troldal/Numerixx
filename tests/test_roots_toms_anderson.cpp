// Toms748 and AndersonBjorck in nxx::roots (Phase 4, docs/roots/PHASE4_API.md).
//
// Evaluation counts are asserted where no libm is involved (+ - * / and sqrt in double), against an independent Python
// implementation of the specification (.scratch/roots4/proto4.py, revision 2), assuming no FMA contraction in fn or in
// the steps, as in Phase 2. The paper's test problems are written with their powers as products, as in the prototype.
// Every other case checks accuracy (roots_check.hpp).

#include "roots_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>

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

    // Wilkinson's W10 = (x - 1) ... (x - 10), expanded, by Horner's scheme: a root at 7 whose computed values are
    // rounding noise within about 1e-10 of it.
    constexpr auto w10 = [](double x) {
        constexpr double c[] = { 1.0, -55.0, 1320.0, -18150.0, 157773.0, -902055.0, 3416930.0, -8409500.0, 12753576.0,
            -10628640.0, 3628800.0 };
        double p = 0.0;
        for (double const a : c) p = p * x + a;
        return p;
    };

    constexpr double cycle_root = -1.769292354238631415240409; // x^3 - 2x + 2 (mpmath)
    constexpr double kepler_root = 1.862086686874532271834067; // E - 0.9 sin E = 1, 0.9 as a double (mpmath)

    // y^n by repeated multiplication, left to right, as the prototype.
    constexpr double power(double y, int n)
    {
        double r = y;
        for (int i = 1; i < n; ++i) r = r * y;
        return r;
    }

    // Alefeld, Potra and Shi (1995), Table I.
    auto problem2(int n)
    {
        return [](double x) {
            double s = 0.0;
            for (int i = 1; i <= 20; ++i) {
                double const t = x - double(i * i);
                s = s + double((2 * i - 5) * (2 * i - 5)) / (t * t * t);
            }
            return -2.0 * s;
        };
        static_cast<void>(n);
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
    auto problem11(int n)
    {
        return [n](double x) { return (double(n) * x - 1.0) / (double(n - 1) * x); };
    }

    // The bound's meaning: fn changes sign within abs_error of value (W10's check in test_roots_open.cpp).
    template<typename F>
    bool sign_change_within(F const &fn, rt::Root<double> const &r)
    {
        double const at = fn(r.value);
        double const lo = fn(r.value - r.abs_error);
        double const hi = fn(r.value + r.abs_error);
        return at == 0.0 || lo == 0.0 || hi == 0.0 || (lo < 0.0) != (at < 0.0) || (hi < 0.0) != (at < 0.0);
    }

    using B = rt::Bracket<double>;

} // namespace

// Constant evaluation, per method and type, and the wide brackets in double.
static_assert(nxx::math::abs(rt::solve<rt::Toms748>(sq, B{ 1.0, 2.0 })->value - sqrt2) <= 3e-16);
static_assert(nxx::math::abs(rt::solve<rt::AndersonBjorck>(sq, B{ 1.0, 2.0 })->value - sqrt2) <= 3e-16);
static_assert(nxx::math::abs(rt::solve<rt::Toms748>([](float x) { return x * x - 2.0f; }, rt::Bracket{ 1.0f, 2.0f })
                                 ->value
                             - 1.41421356f)
              <= 4e-7f);
static_assert(
    nxx::math::abs(rt::solve<rt::AndersonBjorck>([](float x) { return x * x - 2.0f; }, rt::Bracket{ 1.0f, 2.0f })
                       ->value
                   - 1.41421356f)
    <= 4e-7f);
static_assert(
    nxx::math::abs(rt::solve<rt::Toms748>([](long double x) { return x * x - 2.0L; }, rt::Bracket{ 1.0L, 2.0L })->value
                   - 1.4142135623730950488016887242096980785696718753769L)
    <= 8 * std::numeric_limits<long double>::epsilon()); // 2 tol + c, about 7 eps near sqrt(2)
static_assert(nxx::math::abs(rt::solve<rt::AndersonBjorck>([](long double x) { return x * x - 2.0L; },
                                 rt::Bracket{ 1.0L, 2.0L })
                                 ->value
                             - 1.4142135623730950488016887242096980785696718753769L)
              <= 8 * std::numeric_limits<long double>::epsilon()); // 2 tol + c, about 7 eps near sqrt(2)
static_assert(rt::solve<rt::Toms748>(minus_1, B{ -dmax, dmax })->value == 1.0);
static_assert(rt::solve<rt::AndersonBjorck>(minus_1, B{ -dmax, dmax })->value == 1.0);

TEST_CASE("roots: Toms748 and AndersonBjorck, evaluation counts in double")
{
    auto const check = [](auto const &fn, B const &bracket, int toms, int ab) {
        CHECK(counted<rt::Toms748>(fn, bracket).second == toms);
        CHECK(counted<rt::AndersonBjorck>(fn, bracket).second == ab);
    };
    check(sq, { 1.0, 2.0 }, 9, 9);
    check(wal, { 2.0, 3.0 }, 9, 9);
    check(x20, { 0.0, 5.0 }, 21, 30); // AndersonBjorck's guard: 65,109,847 evaluations without it
    check(triple, { 0.0, 3.0 }, 149, 152);
    check(tiny_root, { 0.0, 1.0 }, 3, 4);
    check(minus_1, { -dmax, dmax }, 8, 13); // the wide guard, after zero first
    check(minus_1, { 0.0, dmax }, 7, 12);
    check(cycle, { -2.0, -1.0 }, 9, 10);
    check(reciprocal, { 0.5, 3.0 }, 6, 5);
    check(w10, { 6.7, 7.6 }, 18, 18);
    check(jump, { -1.0, 2.0 }, 981, 51); // the walk to 0: AndersonBjorck's is geometric in the exponent
    check(sqrt_fn, { 0.0, 1.0 }, 6, 9);
    // The pole: AndersonBjorck returns it as a root; Toms748 lands on fl(1/3), where fn is infinite.
    auto const [toms_pole, n_toms_pole] = counted<rt::Toms748>(pole, B{ 0.0, 1.0 });
    CHECK(toms_pole == std::unexpected(rt::error::non_finite));
    CHECK(n_toms_pole == 84);
    CHECK(counted<rt::AndersonBjorck>(pole, B{ 0.0, 1.0 }).second == 95);

    // Alefeld, Potra and Shi's problems.
    for (auto const &[n, toms, ab] : { std::tuple{ 1, 16, 23 }, std::tuple{ 5, 13, 17 }, std::tuple{ 10, 13, 18 } })
        check(problem2(n), { n * n + 1e-9, (n + 1) * (n + 1) - 1e-9 }, toms, ab);
    for (auto const &[n, toms, ab] : { std::tuple{ 4, 13, 19 }, std::tuple{ 8, 16, 22 }, std::tuple{ 12, 18, 23 } })
        check(problem4(n), { 0.0, 5.0 }, toms, ab);
    for (auto const &[n, toms, ab] : { std::tuple{ 5, 9, 9 }, std::tuple{ 10, 7, 8 }, std::tuple{ 20, 8, 8 } })
        check(problem7(n), { 0.0, 1.0 }, toms, ab);
    for (auto const &[n, toms, ab] : { std::tuple{ 5, 11, 10 }, std::tuple{ 10, 11, 11 }, std::tuple{ 20, 12, 13 } })
        check(problem8(n), { 0.0, 1.0 }, toms, ab);
    for (auto const &[n, toms, ab] : { std::tuple{ 5, 11, 10 }, std::tuple{ 15, 9, 8 }, std::tuple{ 20, 9, 8 } })
        check(problem9(n), { 0.0, 1.0 }, toms, ab);
    for (auto const &[n, toms, ab] : { std::tuple{ 5, 19, 7 }, std::tuple{ 15, 20, 6 }, std::tuple{ 20, 18, 6 } })
        check(problem11(n), { 0.01, 1.0 }, toms, ab);
}

TEST_CASE("roots: Toms748 and AndersonBjorck, accuracy in double")
{
    auto const both = [](auto const &fn, B const &bracket, double ref, double c) {
        CHECK(brackets(rt::solve<rt::Toms748>(fn, bracket), ref, c));
        CHECK(brackets(rt::solve<rt::AndersonBjorck>(fn, bracket), ref, c));
    };
    // The conditioning terms of test_roots_bracketing.cpp and test_roots_open.cpp: x^2 - 2 7.9e-17; Wallis 2e-16;
    // x^20 - 1 near 1 1.1e-16; x - 1, (x - 1)^3 and x - 1e-300 exact; x^3 - 2x + 2 1.8e-16; 1/x - 1 1.2e-16; sqrt
    // correctly rounded, exactly 0 at 0.25.
    both(sq, { 1.0, 2.0 }, sqrt2, 8e-17);
    both(wal, { 2.0, 3.0 }, wallis, 2e-16);
    both(x20, { 0.0, 5.0 }, 1.0, 1.1e-16);
    both(triple, { 0.0, 3.0 }, 1.0, 0.0);
    both(tiny_root, { 0.0, 1.0 }, 1e-300, 0.0);
    both(minus_1, { -dmax, dmax }, 1.0, 0.0);
    both(cycle, { -2.0, -1.0 }, cycle_root, 1.8e-16);
    both(reciprocal, { 0.5, 3.0 }, 1.0, 1.2e-16);
    both(sqrt_fn, { 0.0, 1.0 }, 0.25, 0.0);
    // The jump at 0: the sign change is at 0 exactly.
    for (auto const &r : { rt::solve<rt::Toms748>(jump, B{ -1.0, 2.0 }), rt::solve<rt::AndersonBjorck>(jump, B{ -1.0, 2.0 }) }) {
        REQUIRE(r.has_value());
        CHECK(nxx::math::abs(r->value) <= r->abs_error);
    }
    // libm, with Phase 1's budget of 2 ulps: cos x - x 1.4e-16; e^x - 2 4.5e-16; x e^x - 1 1.3e-16; Kepler 2.9e-16;
    // ln(x / 1e-12) 2.5e-28.
    both([](double x) { return std::cos(x) - x; }, { 0.0, 1.0 }, dottie, 1.4e-16);
    both([](double x) { return std::exp(x) - 2.0; }, { 0.0, 1.0 }, ln2, 4.5e-16);
    both([](double x) { return x * std::exp(x) - 1.0; }, { 0.0, 1.0 }, omega, 1.3e-16);
    both([](double x) { return x - 0.9 * std::sin(x) - 1.0; }, { 0.0, 3.141592653589793 }, kepler_root, 2.9e-16);
    both([](double x) { return std::log(x / 1e-12); }, { 1e-20, 1.0 }, 1e-12, 2.5e-28);

    // W10 and the paper's problems, whose computed values are noisy or whose roots have no closed form: the bound's
    // meaning, a sign change within abs_error of value.
    for (auto const &r : { rt::solve<rt::Toms748>(w10, B{ 6.7, 7.6 }), rt::solve<rt::AndersonBjorck>(w10, B{ 6.7, 7.6 }) }) {
        REQUIRE(r.has_value());
        CHECK(sign_change_within(w10, *r));
    }
    auto const sign_change = [](auto const &fn, B const &bracket) {
        for (auto const &r : { rt::solve<rt::Toms748>(fn, bracket), rt::solve<rt::AndersonBjorck>(fn, bracket) }) {
            REQUIRE(r.has_value());
            CHECK(sign_change_within(fn, *r));
            CHECK(r->abs_error <= 2.0 * default_tol(r->value));
        }
    };
    for (int const n : { 1, 5, 10 }) sign_change(problem2(n), { n * n + 1e-9, (n + 1) * (n + 1) - 1e-9 });
    for (int const n : { 5, 10, 20 }) {
        sign_change(problem7(n), { 0.0, 1.0 });
        sign_change(problem8(n), { 0.0, 1.0 });
    }
    for (int const n : { 5, 15, 20 }) {
        sign_change(problem9(n), { 0.0, 1.0 });
        sign_change(problem11(n), { 0.01, 1.0 });
    }
}

TEST_CASE("roots: Toms748 and AndersonBjorck, edges")
{
    // fn exactly 0 at an interior point: the first secant lands on 0.5.
    for (auto const &r : { rt::solve<rt::Toms748>([](double x) { return x - 0.5; }, B{ 0.0, 1.0 }),
             rt::solve<rt::AndersonBjorck>([](double x) { return x - 0.5; }, B{ 0.0, 1.0 }) })
        CHECK(r == rt::Root<double>{ 0.5, default_tol(0.5) });

    // A reversed bracket gives the same result; a tuple fn gives the value fn's, bit for bit.
    auto const pair_fn = [](double x) { return std::pair{ (x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0 }; };
    CHECK(rt::solve<rt::Toms748>(wal, B{ 3.0, 2.0 }) == rt::solve<rt::Toms748>(wal, B{ 2.0, 3.0 }));
    CHECK(rt::solve<rt::AndersonBjorck>(wal, B{ 3.0, 2.0 }) == rt::solve<rt::AndersonBjorck>(wal, B{ 2.0, 3.0 }));
    CHECK(rt::solve<rt::Toms748>(pair_fn, B{ 2.0, 3.0 }) == rt::solve<rt::Toms748>(wal, B{ 2.0, 3.0 }));
    CHECK(rt::solve<rt::AndersonBjorck>(pair_fn, B{ 2.0, 3.0 }) == rt::solve<rt::AndersonBjorck>(wal, B{ 2.0, 3.0 }));

    // The placement's two fallbacks.
    // 1. 2 delta rounds away at an end of much larger magnitude than u: u = 0.3, tol(0.3) = 1.3e-16, and
    //    fl(2.5 - 2 delta) = 2.5; the point becomes the midpoint.
    rt::detail::Tolerance<double> const defaults{ dmin, 2.0 * deps };
    rt::detail::Enclosure<double> const far{ 0.3, -1e-20, 2.5, 1.0, 2.5, 1.0, 2.5, 1.0 };
    CHECK(rt::detail::place(far, std::optional<double>(2.5 - 1e-16), defaults) == rt::detail::midpoint(0.3, 2.5));
    // 2. delta flushes to 0 where denorm_min() is min(), as in cpp_bin_float_50: at u = 0 delta is min, at least.
    rt::detail::Tolerance<bin50> const big_defaults{ bin50(std::numeric_limits<bin50>::denorm_min()),
        bin50(2 * std::numeric_limits<bin50>::epsilon()) };
    bin50 const smallest = std::numeric_limits<bin50>::denorm_min();
    rt::detail::Enclosure<bin50> const near0{ bin50(0), bin50(-smallest), bin50(0.5), bin50(1), bin50(0.5), bin50(1),
        bin50(0.5), bin50(1) };
    bin50 const placed = rt::detail::place(near0, std::optional<bin50>(bin50(0)), big_defaults);
    CHECK(placed > bin50(0));
    CHECK(placed < bin50(0.5));
    // ... and the solve that reached it: fn = 1 above 0 and -min at 0, on [0, 0.5].
    auto const step_at_0 = [smallest](bin50 const &x) { return x > 0 ? bin50(1) : bin50(-smallest); };
    for (auto const &r : { rt::solve<rt::Toms748>(step_at_0, rt::Bracket<bin50>{ bin50(0), bin50(0.5) }),
             rt::solve<rt::AndersonBjorck>(step_at_0, rt::Bracket<bin50>{ bin50(0), bin50(0.5) }) }) {
        REQUIRE(r.has_value());
        CHECK(r->abs_error <= bin50(4 * smallest));
    }
}

TEST_CASE("roots: Toms748 and AndersonBjorck in float, long double and multiprecision types")
{
    // x^2 - 2 and the Wallis cubic with c = eps |ref| (fn's rounding), as for Brent.
    auto const check = [](auto zero) {
        using T = decltype(zero);
        T const eps = std::numeric_limits<T>::epsilon();
        T const tiny = T(std::numeric_limits<T>::denorm_min());
        T const s2 = from_string<T>(sqrt2_70);
        T const w = from_string<T>(wallis_70);
        auto const sq_t = [](T const &x) { return T(x * x - 2); };
        auto const cubic = [](T const &x) { return T(T(x * x - 2) * x - 5); };
        for (auto const &r : { rt::solve<rt::Toms748>(sq_t, rt::Bracket<T>{ T(1), T(2) }),
                 rt::solve<rt::AndersonBjorck>(sq_t, rt::Bracket<T>{ T(1), T(2) }) })
            CHECK(brackets(r, s2, T(eps * s2), T(2 * T(tiny + T(2 * eps * s2)))));
        for (auto const &r : { rt::solve<rt::Toms748>(cubic, rt::Bracket<T>{ T(2), T(3) }),
                 rt::solve<rt::AndersonBjorck>(cubic, rt::Bracket<T>{ T(2), T(3) }) })
            CHECK(brackets(r, w, T(eps * w), T(2 * T(tiny + T(2 * eps * w)))));
    };
    check(float{});
    check(static_cast<long double>(0));
    check(bin50{});
    check(bin50_et{});
    check(dec50{});
}
