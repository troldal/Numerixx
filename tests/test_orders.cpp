// Derivatives of order 1 to 4 (Phase 4): the order argument, the generated tables, exactness on polynomials, overflow
// before the final quotient, and Ridders for f'' to f''''. References are closed forms, published tables (Fornberg
// 1988, Chapra ch. 23) and the moment conditions computed here; Ridders' per-case values come from the independent
// Python prototype (docs/PHASE4_API.md), which the C++ results for double match exactly.

#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <numeric>
#include <type_traits>

namespace mp = boost::multiprecision;

using bin50 = mp::cpp_bin_float_50;
using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
using dec50 = mp::cpp_dec_float_50;

using nxx::deriv::Backward;
using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Estimate;
using nxx::deriv::Forward;
using nxx::deriv::order;
using nxx::deriv::Ridders;

namespace {

    template<typename S, int N>
    constexpr auto table = nxx::deriv::detail::make_table<N>(S{});

    // M_j = sum_k numerators[k] n_k^j, computed here: a Taylor check, independent of the interpolation formula of the
    // generator. The largest term is about 6.2e6 * 6^14 = 4.8e17 (Central<10>, N = 4, j = 14), inside long long.
    template<std::size_t Size>
    constexpr long long moment(nxx::deriv::detail::Table<Size> const &t, int j)
    {
        long long m = 0;
        for (std::size_t k = 0; k < Size; ++k) {
            long long power = 1;
            for (int i = 0; i < j; ++i) power *= t.nodes[k];
            m += t.numerators[k] * power;
        }
        return m;
    }

    constexpr long long factorial(int n) { return n <= 1 ? 1 : n * factorial(n - 1); }

    // A rule of accuracy order p for the N-th derivative has M_j = 0 for j != N up to N + p - 1, M_N = N! den and
    // M_(N+p) != 0 (so the order is exactly p), strictly increasing nodes, a reduced denominator (with the moments
    // that fixes the table uniquely) and values below 2^24, exact in float.
    template<typename S, int N>
    constexpr bool valid_table()
    {
        constexpr auto t = table<S, N>;
        constexpr int p = t.order;
        for (std::size_t k = 1; k < t.nodes.size(); ++k)
            if (!(t.nodes[k - 1] < t.nodes[k])) return false;
        for (int j = 0; j < N + p; ++j)
            if (moment(t, j) != (j == N ? factorial(N) * t.denominator : 0)) return false;
        if (moment(t, N + p) == 0) return false;
        int g = t.denominator;
        for (int const v : t.numerators) g = std::gcd(g, v);
        if (g != 1 || t.denominator <= 0 || t.denominator >= (1 << 24)) return false;
        for (int const v : t.numerators)
            if (v >= (1 << 24) || v <= -(1 << 24)) return false;
        return true;
    }

    template<int N>
    constexpr bool valid_tables()
    {
        return valid_table<Central<2>, N>() && valid_table<Central<4>, N>() && valid_table<Central<6>, N>()
               && valid_table<Central<8>, N>() && valid_table<Central<10>, N>() && valid_table<Forward<1>, N>()
               && valid_table<Forward<2>, N>() && valid_table<Forward<3>, N>() && valid_table<Forward<4>, N>()
               && valid_table<Forward<5>, N>() && valid_table<Forward<6>, N>() && valid_table<Backward<1>, N>()
               && valid_table<Backward<2>, N>() && valid_table<Backward<3>, N>() && valid_table<Backward<4>, N>()
               && valid_table<Backward<5>, N>() && valid_table<Backward<6>, N>();
    }

    // Backward<P> is Forward<P> with h replaced by -h: the nodes mirrored, the weights times (-1)^N.
    template<int P, int N>
    constexpr bool mirrored()
    {
        constexpr auto f = table<Forward<P>, N>;
        constexpr auto b = table<Backward<P>, N>;
        constexpr std::size_t size = f.nodes.size();
        constexpr int sign = N % 2 == 1 ? -1 : 1;
        for (std::size_t k = 0; k < size; ++k)
            if (b.nodes[k] != -f.nodes[size - 1 - k] || b.numerators[k] != sign * f.numerators[size - 1 - k])
                return false;
        return b.denominator == f.denominator;
    }

    template<typename S, int N, std::size_t Size>
    constexpr bool has_weights(std::array<int, Size> const &nodes, std::array<int, Size> const &numerators, int den)
    {
        constexpr auto t = table<S, N>;
        if constexpr (t.nodes.size() != Size)
            return false;
        else
            return t.nodes == nodes && t.numerators == numerators && t.denominator == den;
    }

    template<int D>
    constexpr auto power = [](double x) {
        double y = 1;
        for (int i = 0; i < D; ++i) y *= x;
        return y;
    };

    // The N-th derivative of x^D at 1: D (D - 1) ... (D - N + 1).
    constexpr double falling(int d, int n)
    {
        double r = 1;
        for (int i = 0; i < n; ++i) r *= d - i;
        return r;
    }

    // Whether x^D under the rule stays exact in double at x = a / 4, h = b / 4: every value (a + n_k b)^D / 4^D,
    // every term and every partial sum must fit 53 bits, which sum_k |numerators[k]| |a + n_k b|^D < 2^53 ensures
    // (computed in double, which only rounds once the sum is far above the bound).
    template<typename S, int N, int D>
    constexpr bool fits_53_bits(int a, int b)
    {
        constexpr auto t = table<S, N>;
        double bound = 0;
        for (std::size_t k = 0; k < t.nodes.size(); ++k) {
            int const base = a + t.nodes[k] * b;
            double value = 1;
            for (int i = 0; i < D; ++i) value *= base < 0 ? -base : base;
            bound += (t.numerators[k] < 0 ? -t.numerators[k] : t.numerators[k]) * value;
        }
        return bound < 9007199254740992.0;
    }

    // A rule of order p for the N-th derivative is exact for x^(N + p - 1) at x = A/4 with h = 1/4, when the bit
    // budget allows: A = 4 (x = 1) for every rule but Central<10> for N = 3 and 4, which need A = 2.
    template<typename S, int N, int A = 4>
    constexpr bool exact_on_polynomial()
    {
        constexpr int d = N + table<S, N>.order - 1;
        static_assert(fits_53_bits<S, N, d>(A, 1));
        double const x = A / 4.0;
        double expected = falling(d, N);
        for (int i = 0; i < d - N; ++i) expected *= x;
        return derivative<S, N>(power<d>, x, 0.25)->value == expected;
    }

    template<int N>
    constexpr bool all_exact()
    {
        constexpr int wide = N >= 3 ? 2 : 4;
        return exact_on_polynomial<Central<2>, N>() && exact_on_polynomial<Central<4>, N>()
               && exact_on_polynomial<Central<6>, N>() && exact_on_polynomial<Central<8>, N>()
               && exact_on_polynomial<Central<10>, N, wide>() && exact_on_polynomial<Forward<1>, N>()
               && exact_on_polynomial<Forward<2>, N>() && exact_on_polynomial<Forward<3>, N>()
               && exact_on_polynomial<Forward<4>, N>() && exact_on_polynomial<Forward<5>, N>()
               && exact_on_polynomial<Forward<6>, N>() && exact_on_polynomial<Backward<1>, N>()
               && exact_on_polynomial<Backward<2>, N>() && exact_on_polynomial<Backward<3>, N>()
               && exact_on_polynomial<Backward<4>, N>() && exact_on_polynomial<Backward<5>, N>()
               && exact_on_polynomial<Backward<6>, N>();
    }

    constexpr double dmax = std::numeric_limits<double>::max();
    constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
    constexpr auto non_finite = std::unexpected(nxx::deriv::error::non_finite);

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

} // namespace

// The order: an integer or an order, 1 to 4, and nothing else.
template<typename M, auto N>
concept accepted = requires { derivative<M, N>(power<2>, 1.0); };
static_assert(accepted<Central<2>, 1> && accepted<Central<2>, 4> && accepted<Ridders, order::third>);
static_assert(accepted<Forward<6>, 4u> && accepted<Backward<1>, std::size_t{ 2 }> && accepted<Central<10>, 3L>);
static_assert(!accepted<Central<2>, 0> && !accepted<Central<2>, 5> && !accepted<Ridders, -1>);
static_assert(!accepted<Central<2>, true> && !accepted<Central<2>, 2.0>);
static_assert(!accepted<Central<2>, 0x1'0000'0001ULL>); // not wrapped into order 1
// The low-precision warning: a 4th derivative in float only.
static_assert(nxx::deriv::detail::noisy<4, float> && nxx::deriv::detail::noisy<order::fourth, float>);
static_assert(!nxx::deriv::detail::noisy<4, double> && !nxx::deriv::detail::noisy<3, float>);
static_assert(!nxx::deriv::detail::noisy<4, long double> && !nxx::deriv::detail::noisy<4, bin50>);
static_assert(static_cast<int>(order::first) == 1 && static_cast<int>(order::fourth) == 4);
static_assert(!std::is_convertible_v<order, int>); // scoped: no implicit conversion elsewhere
static_assert(
    derivative<Central<2>, order::second>(power<3>, 1.0, 0.25) == derivative<Central<2>, 2>(power<3>, 1.0, 0.25));
static_assert(derivative<Central<4>, order::first>(power<3>, 1.0, 0.25) == derivative<Central<4>>(power<3>, 1.0, 0.25));

// Every table of every order satisfies its moment conditions and fits 24 bits.
static_assert(valid_tables<1>() && valid_tables<2>() && valid_tables<3>() && valid_tables<4>());
static_assert(mirrored<1, 1>() && mirrored<2, 2>() && mirrored<3, 3>() && mirrored<4, 4>() && mirrored<6, 4>());
static_assert(mirrored<6, 1>() && mirrored<5, 2>() && mirrored<1, 3>() && mirrored<2, 4>());

// Published rows: Fornberg (1988), Tables 1 (centred) and 3 (one-sided), and Chapra ch. 23 (the backward f'''' of
// O(h^2), 3 f(x) - 14 f(x - h) + 26 f(x - 2h) - 24 f(x - 3h) + 11 f(x - 4h) - 2 f(x - 5h)).
static_assert(has_weights<Central<2>, 2, 3>({ -1, 0, 1 }, { 1, -2, 1 }, 1));
static_assert(has_weights<Central<4>, 2, 5>({ -2, -1, 0, 1, 2 }, { -1, 16, -30, 16, -1 }, 12));
static_assert(has_weights<Central<6>, 2, 7>({ -3, -2, -1, 0, 1, 2, 3 }, { 2, -27, 270, -490, 270, -27, 2 }, 180));
static_assert(has_weights<Central<2>, 3, 4>({ -2, -1, 1, 2 }, { -1, 2, -2, 1 }, 2));
static_assert(has_weights<Central<4>, 3, 6>({ -3, -2, -1, 1, 2, 3 }, { 1, -8, 13, -13, 8, -1 }, 8));
static_assert(has_weights<Central<2>, 4, 5>({ -2, -1, 0, 1, 2 }, { 1, -4, 6, -4, 1 }, 1));
static_assert(has_weights<Central<4>, 4, 7>({ -3, -2, -1, 0, 1, 2, 3 }, { -1, 12, -39, 56, -39, 12, -1 }, 6));
static_assert(has_weights<Forward<2>, 2, 4>({ 0, 1, 2, 3 }, { 2, -5, 4, -1 }, 1));
static_assert(has_weights<Forward<3>, 2, 5>({ 0, 1, 2, 3, 4 }, { 35, -104, 114, -56, 11 }, 12));
static_assert(has_weights<Forward<1>, 3, 4>({ 0, 1, 2, 3 }, { -1, 3, -3, 1 }, 1));
static_assert(has_weights<Forward<2>, 3, 5>({ 0, 1, 2, 3, 4 }, { -5, 18, -24, 14, -3 }, 2));
static_assert(has_weights<Forward<1>, 4, 5>({ 0, 1, 2, 3, 4 }, { 1, -4, 6, -4, 1 }, 1));
static_assert(has_weights<Backward<2>, 4, 6>({ -5, -4, -3, -2, -1, 0 }, { -2, 11, -24, 26, -14, 3 }, 1));

// Exactness on polynomials of the highest degree a rule integrates, x^(N + p - 1), for every rule and order.
static_assert(all_exact<1>() && all_exact<2>() && all_exact<3>() && all_exact<4>());
// float, for the low orders, whose values fit 24 bits: x^3 under Central<2> for f'' and x^2 under Forward<1> for f''
// at 1 with h = 1/4 (at most 5^3 = 125 in units of 4^-3).
static_assert(derivative<Central<2>, 2>([](float x) { return x * x * x; }, 1.0f, 0.25f)->value == 6.0f);
static_assert(derivative<Forward<1>, 2>([](float x) { return x * x; }, 1.0f, 0.25f)->value == 2.0f);

// The actual offsets of Central<2> for f'': at x = 1 with h = 0x1.4p-53 the points round to 1 + 2^-52 and 1 - 2^-53,
// unequal offsets. For fn(t) = t, f+ - f0 and f0 - f- are exactly the offsets, so f'' = 0 exactly; a formula with the
// nominal h would give 2^-53 / h^2 = 5.8e15.
static_assert(derivative<Central<2>, 2>(power<1>, 1.0, 0x1.4p-53) == 0.0);

// h^N is never formed: each step below is outside the range of h^N, and the results are exact. N = 2: f(+-h) =
// 2^-600, so (1 - -1) / 2^-600 = 2^601.
namespace {
    constexpr auto wide2 = [](double t) { return (t * 0x1p600) * t; };
    constexpr auto narrow2 = [](double t) { return (t * 0x1p-600) * t; };
    constexpr auto wide3 = [](double t) { return (0x1p300 * t) * (0x1p300 * t) * (0x1p300 * t); };
    constexpr auto narrow3 = [](double t) { return (t * 0x1p-260) * (t * 0x1p-260) * (t * 0x1p-260); };
    constexpr auto wide4 = [](double t) { return power<4>(0x1p100 * t); };
    constexpr auto narrow4 = [](double t) { return power<4>(t * 0x1p-100); };
} // namespace
static_assert(derivative<Central<2>, 2>(wide2, 0.0, 0x1p-600) == 0x1p601);
static_assert(derivative<Central<2>, 2>(narrow2, 0.0, 0x1p600) == 0x1p-599);
static_assert(derivative<Central<2>, 3>(wide3, 0.0, 0x1p-600) == 6 * 0x1p900);
static_assert(derivative<Central<2>, 3>(narrow3, 0.0, 0x1p600) == 6 * 0x1p-780);
static_assert(derivative<Central<2>, 4>(wide4, 0.0, 0x1p-300) == 24 * 0x1p400);
static_assert(derivative<Central<2>, 4>(narrow4, 0.0, 0x1p300) == 24 * 0x1p-400);

// Overflow before the final quotient is an error value, also in constant expressions. Regression: a constant 1e308
// under Central<4> (f') formed -8e308 and then inf - inf, which GCC and Clang rejected at compile time; the true
// derivative is 0, but the weighted sum cannot be formed.
namespace {
    constexpr auto big = [](double) { return 1e308; };
    // f(0) = -max, f(+-h) = max: f+ - f0 = 2 max.
    constexpr auto spike = [](double t) { return t == 0.0 ? -dmax : dmax; };
    // f(+h) = 1e300, f(0) = 0, f(-h) = -1e300: with h = 1e-10 both one-sided quotients overflow, with the same sign.
    constexpr auto jump = [](double t) { return t > 0.0 ? 1e300 : (t < 0.0 ? -1e300 : 0.0); };
} // namespace
static_assert(derivative<Central<4>>(big, 0.0, 1.0) == non_finite);
static_assert(derivative<Central<2>, 4>(big, 0.0, 1.0) == non_finite);
static_assert(derivative<Forward<3>, 3>(big, 0.0, 1.0) == non_finite);
static_assert(derivative<Central<2>, 2>(spike, 0.0, 1.0) == non_finite);
static_assert(derivative<Central<2>, 2>(jump, 0.0, 1e-10) == non_finite);
// Regression: fl(max / 6) rounds up, and 6 times it ties to infinity, so the term check must be >=, not >.
namespace {
    constexpr auto sixth = [](double) { return dmax / 6; };
} // namespace
static_assert(derivative<Central<2>, 4>(sixth, 0.0, 1.0) == non_finite);
static_assert(derivative<Ridders, 4>(sixth, 0.0, 1.0) == non_finite);
static_assert(derivative<Ridders, 2>(jump, 0.0, 1e-10) == non_finite);

// Ridders on polynomials, x = 2 and h = 1/8, every value dyadic: x^(N+1) is exact at level 1, where the base is exact
// (5, 8 and 9 evaluations); x^(N+2) needs column 1 of the extrapolation (7, 12 and 13 evaluations).
static_assert(derivative<Ridders, 2>(power<3>, 2.0, 0.125)->value == 12.0);
static_assert(derivative<Ridders, 3>(power<4>, 2.0, 0.125)->value == 48.0);
static_assert(derivative<Ridders, 4>(power<5>, 2.0, 0.125)->value == 240.0);
static_assert(derivative<Ridders, 2>(power<4>, 2.0, 0.125)->value == 48.0);
static_assert(derivative<Ridders, 3>(power<5>, 2.0, 0.125)->value == 240.0);
static_assert(derivative<Ridders, 4>(power<6>, 2.0, 0.125)->value == 1440.0);

// The floor saturates at max at every step: a constant has every base value exactly 0, so the tableau stops on exact
// agreement with the floor as the estimate. fn = 1e300 with h = 1e-300 overflows the first division of the floor,
// fn = 1 only the last (h = 1e-170, 1e-130 and 1e-100 for N = 2, 3, 4).
namespace {
    constexpr auto huge = [](double) { return 1e300; };
    constexpr auto unit = [](double) { return 1.0; };
    constexpr Estimate<double> saturated{ .value = 0.0, .abs_error = dmax };
} // namespace
static_assert(derivative<Ridders, 2>(huge, 0.0, 1e-300) == saturated);
static_assert(derivative<Ridders, 3>(huge, 0.0, 1e-300) == saturated);
static_assert(derivative<Ridders, 4>(huge, 0.0, 1e-300) == saturated);
static_assert(derivative<Ridders, 2>(unit, 0.0, 1e-170) == saturated);
static_assert(derivative<Ridders, 3>(unit, 0.0, 1e-130) == saturated);
static_assert(derivative<Ridders, 4>(unit, 0.0, 1e-100) == saturated);

// fn(x) not finite: every rule with a node at x fails, odd-order central rules and Ridders never evaluate x.
namespace {
    constexpr auto hole = [](double t) { return t == 0.0 ? qnan : 1.0; };
} // namespace
static_assert(derivative<Central<2>, 2>(hole, 0.0, 0.5) == non_finite);
static_assert(derivative<Ridders, 2>(hole, 0.0) == non_finite);
static_assert(derivative<Ridders, 4>(hole, 0.0) == non_finite);
static_assert(derivative<Central<2>, 3>(hole, 0.0, 0.5) == 0.0);
static_assert(derivative<Ridders, 3>(hole, 0.0, 0.5)->value == 0.0);

// The default steps: eps^(1/(p+N)) max(|x|, 1) for a stencil, max(|x|, 1) / 4 for Ridders of order 2 to 4.
static_assert(derivative<Central<2>, 2>(power<5>, 3.0)
              == derivative<Central<2>, 2>(power<5>, 3.0, nxx::math::root_epsilon<4, double>() * 3.0));
static_assert(derivative<Forward<2>, 3>(power<5>, 0.5)
              == derivative<Forward<2>, 3>(power<5>, 0.5, nxx::math::root_epsilon<5, double>()));
static_assert(derivative<Central<4>, 4>(power<7>, -2.0)
              == derivative<Central<4>, 4>(power<7>, -2.0, nxx::math::root_epsilon<8, double>() * 2.0));
static_assert(derivative<Ridders, 2>(power<5>, 3.0) == derivative<Ridders, 2>(power<5>, 3.0, 0.75));
static_assert(derivative<Ridders, 3>(power<5>, 3.0) == derivative<Ridders, 3>(power<5>, 3.0, 0.75));
static_assert(derivative<Ridders, 4>(power<6>, 0.5) == derivative<Ridders, 4>(power<6>, 0.5, 0.25));
static_assert(derivative<Ridders>(power<5>, 3.0) == derivative<Ridders>(power<5>, 3.0, 3.0 / 32));

namespace {

    // The prototype's values for a case, relative to |f^(N)|, and its number of evaluations.
    struct Reference
    {
        double error;
        double estimate;
        int evaluations;
    };

    // As in test_ridders.cpp: the actual error within `cover` times the estimate (6, from the prototype's worst case
    // of 3.6 over 5,400 random cases for f'; for f'' to f'''' the worst was 2.0), within 4 times the prototype's
    // error, the estimate within twice the prototype's, and the evaluations within the prototype's plus 4: libm
    // implementations round fn's values differently, which can move a stall by a level.
    template<int N, typename T, typename F>
    bool accurate(F fn, T x, T exact, Reference ref, T cover = T(6))
    {
        int calls = 0;
        auto const r = derivative<Ridders, N>(Counted{ fn, &calls }, x);
        if (!r || !r->abs_error) return false;
        T const actual = std::abs(r->value - exact);
        bool const ok = actual <= cover * *r->abs_error && actual <= T(4) * T(ref.error) * std::abs(exact)
                        && *r->abs_error <= T(2) * T(ref.estimate) * std::abs(exact) && calls <= ref.evaluations + 4;
        if (!ok) MESSAGE("N ", N, ": actual ", actual, ", estimate ", *r->abs_error, ", calls ", calls);
        return ok;
    }

    // Only coverage, for float cases whose value depends on the libm's rounding of fn.
    template<int N, typename T, typename F>
    bool covered(F fn, T x, T exact, T cover)
    {
        auto const r = derivative<Ridders, N>(fn, x);
        return r && r->abs_error && std::abs(r->value - exact) <= cover * *r->abs_error;
    }

    // Closed forms of the N-th derivative.
    double sin_n(int n, double x) { return std::sin(x + n * std::numbers::pi / 2); }
    double log_n(int n, double x) { return (n % 2 == 1 ? 1 : -1) * double(factorial(n - 1)) / std::pow(x, n); }
    double pole_n(int n, double a) { return -double(factorial(n)) / std::pow(a, n + 1); } // 1/(x - a) at 0
    // 1/(1 + 25 x^2) = Re 1/(1 + 5ix), whose N-th derivative is (-1)^N N! (5i)^N (1 + 5ix)^-(N+1).
    double runge_n(int n, double x)
    {
        std::complex<double> const z(1.0, 5.0 * x);
        std::complex<double> const c = std::pow(std::complex<double>(0.0, 5.0), n) * double(factorial(n));
        return ((n % 2 == 0 ? c : -c) / std::pow(z, n + 1)).real();
    }

} // namespace

TEST_CASE("orders: the compile-time cases at run time")
{
    CHECK(derivative<Central<2>, order::second>(power<3>, 1.0, 0.25) == derivative<Central<2>, 2>(power<3>, 1.0, 0.25));
    CHECK(derivative<Central<4>, order::first>(power<3>, 1.0, 0.25) == derivative<Central<4>>(power<3>, 1.0, 0.25));
    CHECK(all_exact<1>());
    CHECK(all_exact<2>());
    CHECK(all_exact<3>());
    CHECK(all_exact<4>());
    CHECK(derivative<Central<2>, 2>([](float x) { return x * x * x; }, 1.0f, 0.25f)->value == 6.0f);
    CHECK(derivative<Forward<1>, 2>([](float x) { return x * x; }, 1.0f, 0.25f)->value == 2.0f);
    CHECK(derivative<Central<2>, 2>(power<1>, 1.0, 0x1.4p-53) == 0.0);
    CHECK(derivative<Central<2>, 2>(wide2, 0.0, 0x1p-600) == 0x1p601);
    CHECK(derivative<Central<2>, 2>(narrow2, 0.0, 0x1p600) == 0x1p-599);
    CHECK(derivative<Central<2>, 3>(wide3, 0.0, 0x1p-600) == 6 * 0x1p900);
    CHECK(derivative<Central<2>, 3>(narrow3, 0.0, 0x1p600) == 6 * 0x1p-780);
    CHECK(derivative<Central<2>, 4>(wide4, 0.0, 0x1p-300) == 24 * 0x1p400);
    CHECK(derivative<Central<2>, 4>(narrow4, 0.0, 0x1p300) == 24 * 0x1p-400);
    CHECK(derivative<Central<4>>(big, 0.0, 1.0) == non_finite);
    CHECK(derivative<Central<2>, 4>(big, 0.0, 1.0) == non_finite);
    CHECK(derivative<Forward<3>, 3>(big, 0.0, 1.0) == non_finite);
    CHECK(derivative<Central<2>, 2>(spike, 0.0, 1.0) == non_finite);
    CHECK(derivative<Central<2>, 2>(jump, 0.0, 1e-10) == non_finite);
    CHECK(derivative<Ridders, 2>(jump, 0.0, 1e-10) == non_finite);
    CHECK(derivative<Central<2>, 4>(sixth, 0.0, 1.0) == non_finite);
    CHECK(derivative<Ridders, 4>(sixth, 0.0, 1.0) == non_finite);
    CHECK(derivative<Ridders, 2>(power<3>, 2.0, 0.125)->value == 12.0);
    CHECK(derivative<Ridders, 3>(power<4>, 2.0, 0.125)->value == 48.0);
    CHECK(derivative<Ridders, 4>(power<5>, 2.0, 0.125)->value == 240.0);
    CHECK(derivative<Ridders, 2>(power<4>, 2.0, 0.125)->value == 48.0);
    CHECK(derivative<Ridders, 3>(power<5>, 2.0, 0.125)->value == 240.0);
    CHECK(derivative<Ridders, 4>(power<6>, 2.0, 0.125)->value == 1440.0);
    CHECK(derivative<Ridders, 2>(huge, 0.0, 1e-300) == saturated);
    CHECK(derivative<Ridders, 3>(huge, 0.0, 1e-300) == saturated);
    CHECK(derivative<Ridders, 4>(huge, 0.0, 1e-300) == saturated);
    CHECK(derivative<Ridders, 2>(unit, 0.0, 1e-170) == saturated);
    CHECK(derivative<Ridders, 3>(unit, 0.0, 1e-130) == saturated);
    CHECK(derivative<Ridders, 4>(unit, 0.0, 1e-100) == saturated);
    CHECK(derivative<Central<2>, 2>(hole, 0.0, 0.5) == non_finite);
    CHECK(derivative<Ridders, 2>(hole, 0.0) == non_finite);
    CHECK(derivative<Ridders, 4>(hole, 0.0) == non_finite);
    CHECK(derivative<Central<2>, 3>(hole, 0.0, 0.5) == 0.0);
    CHECK(derivative<Ridders, 3>(hole, 0.0, 0.5)->value == 0.0);
    CHECK(derivative<Central<2>, 2>(power<5>, 3.0)
          == derivative<Central<2>, 2>(power<5>, 3.0, nxx::math::root_epsilon<4, double>() * 3.0));
    CHECK(derivative<Forward<2>, 3>(power<5>, 0.5)
          == derivative<Forward<2>, 3>(power<5>, 0.5, nxx::math::root_epsilon<5, double>()));
    CHECK(derivative<Central<4>, 4>(power<7>, -2.0)
          == derivative<Central<4>, 4>(power<7>, -2.0, nxx::math::root_epsilon<8, double>() * 2.0));
    CHECK(derivative<Ridders, 2>(power<5>, 3.0) == derivative<Ridders, 2>(power<5>, 3.0, 0.75));
    CHECK(derivative<Ridders, 3>(power<5>, 3.0) == derivative<Ridders, 3>(power<5>, 3.0, 0.75));
    CHECK(derivative<Ridders, 4>(power<6>, 0.5) == derivative<Ridders, 4>(power<6>, 0.5, 0.25));
    CHECK(derivative<Ridders>(power<5>, 3.0) == derivative<Ridders>(power<5>, 3.0, 3.0 / 32));

    // Coinciding points: at x = 1 with h = 0x1.4p-53 (0.6 ulp above 1), x + h and x + 2h both round to 1 + 2^-52, so
    // the rules that use both are rejected rather than returning a finite wrong value.
    CHECK(derivative<Central<2>, 4>(power<2>, 1.0, 0x1.4p-53) == non_finite);
    CHECK(derivative<Forward<1>, 2>(power<2>, 1.0, 0x1.4p-53) == non_finite);
}

TEST_CASE("orders: Ridders' evaluations")
{
    int calls = 0;
    // The exact stops of the polynomials above.
    CHECK(derivative<Ridders, 2>(Counted{ power<3>, &calls }, 2.0, 0.125)->value == 12.0);
    CHECK(calls == 5);
    calls = 0;
    CHECK(derivative<Ridders, 3>(Counted{ power<4>, &calls }, 2.0, 0.125)->value == 48.0);
    CHECK(calls == 8);
    calls = 0;
    CHECK(derivative<Ridders, 4>(Counted{ power<5>, &calls }, 2.0, 0.125)->value == 240.0);
    CHECK(calls == 9);
    calls = 0;
    CHECK(derivative<Ridders, 2>(Counted{ power<4>, &calls }, 2.0, 0.125)->value == 48.0);
    CHECK(calls == 7);
    calls = 0;
    CHECK(derivative<Ridders, 3>(Counted{ power<5>, &calls }, 2.0, 0.125)->value == 240.0);
    CHECK(calls == 12);
    calls = 0;
    CHECK(derivative<Ridders, 4>(Counted{ power<6>, &calls }, 2.0, 0.125)->value == 1440.0);
    CHECK(calls == 13);

    // fn(x) not finite costs one evaluation for even N; for odd N x is never evaluated.
    calls = 0;
    CHECK(derivative<Ridders, 2>(Counted{ hole, &calls }, 0.0) == non_finite);
    CHECK(calls == 1);
    calls = 0;
    CHECK(derivative<Ridders, 4>(Counted{ hole, &calls }, 0.0) == non_finite);
    CHECK(calls == 1);

    // fn never finite: every level fails at its first evaluation, 32 levels in all.
    auto const nan = [](double) { return std::numeric_limits<double>::quiet_NaN(); };
    calls = 0;
    CHECK(derivative<Ridders, 3>(Counted{ nan, &calls }, 1.0) == non_finite);
    CHECK(calls == 32);

    // x^(N+4) needs column 2. Levels 0 to 2 are exact (the values (16 + k / 2^i)^8 / 8^8 fit 53 bits up to h = 1/32,
    // and the /3 and /15 of the tableau are exact on 26880 + 13440 h^2 + 504 h^4), but from level 3 (h = 1/64) the
    // values round: 129^8 needs 56 bits. So the result carries the rounding of the level-3 base, at most about
    // 1.5 S eps |f| / h^4 = 1.5 * 16 * eps * 2.2^8 * 64^4 = 2.6e-5 (S = 16), and the observed error is 6.9e-7 (2.6e-11
    // relative, as in the prototype). The estimate must cover it.
    auto const r = derivative<Ridders, 4>(power<8>, 2.0, 0.125);
    REQUIRE(r.has_value());
    CHECK(std::abs(r->value - 26880.0) <= 3e-5);
    CHECK(std::abs(r->value - 26880.0) <= *r->abs_error);
    // The same for N = 3, x^7: here the rounding happens to cancel, and the result is exact.
    auto const r3 = derivative<Ridders, 3>(power<7>, 2.0, 0.125);
    REQUIRE(r3.has_value());
    CHECK(std::abs(r3->value - 3360.0) <= 1e-5);
    CHECK(std::abs(r3->value - 3360.0) <= *r3->abs_error);
    auto const r2 = derivative<Ridders, 2>(power<6>, 2.0, 0.125);
    REQUIRE(r2.has_value());
    CHECK(r2->value == 480.0); // exact for N = 2 (the divisions by 3 and 15 happen to be exact here), 9 evaluations
}

TEST_CASE("orders: Ridders' accuracy against closed forms, double and float")
{
    auto const exp = [](double t) { return std::exp(t); };
    auto const sin = [](double t) { return std::sin(t); };
    auto const log = [](double t) { return std::log(t); };
    auto const runge = [](double t) { return 1.0 / (1.0 + 25.0 * t * t); };
    auto const pole = [](double t) { return 1.0 / (t - 0.01); };
    double const e = std::numbers::e;

    // Prototype references {relative error, relative estimate, evaluations}, from the default initial step s/4.
    CHECK(accurate<2>(exp, 1.0, e, { 4.0e-13, 3.6e-12, 15 }));
    CHECK(accurate<2>(sin, 0.5, sin_n(2, 0.5), { 4.1e-13, 3.6e-12, 15 }));
    CHECK(accurate<2>(runge, 0.3, runge_n(2, 0.3), { 2.4e-13, 3.6e-12, 19 }));
    CHECK(accurate<2>(log, 0.2, log_n(2, 0.2), { 5.5e-13, 3.4e-11, 21 }));
    CHECK(accurate<2>(log, 0.01, log_n(2, 0.01), { 9.1e-12, 6.3e-11, 29 })); // domain edge, restarts
    CHECK(accurate<2>(pole, 0.0, pole_n(2, 0.01), { 1.3e-12, 1.2e-11, 55 })); // pole between samples
    CHECK(accurate<3>(exp, 1.0, e, { 3.4e-11, 7.2e-10, 28 }));
    CHECK(accurate<3>(sin, 0.5, sin_n(3, 0.5), { 1.4e-11, 1.4e-9, 28 }));
    CHECK(accurate<3>(runge, 0.3, runge_n(3, 0.3), { 6.9e-12, 6.4e-9, 40 }));
    CHECK(accurate<3>(log, 0.2, log_n(3, 0.2), { 6.0e-11, 4.6e-9, 34 }));
    CHECK(accurate<3>(log, 0.01, log_n(3, 0.01), { 9.9e-10, 6.7e-8, 42 }));
    CHECK(accurate<3>(pole, 0.0, pole_n(3, 0.01), { 2.2e-11, 7.5e-9, 108 }));
    CHECK(accurate<4>(exp, 1.0, e, { 5.4e-10, 2.1e-7, 25 }));
    CHECK(accurate<4>(sin, 0.5, sin_n(4, 0.5), { 3.5e-10, 1.3e-7, 25 }));
    CHECK(accurate<4>(runge, 0.3, runge_n(4, 0.3), { 3.6e-9, 9.9e-7, 37 }));
    CHECK(accurate<4>(log, 0.2, log_n(4, 0.2), { 2.0e-8, 6.0e-6, 35 }));
    CHECK(accurate<4>(log, 0.01, log_n(4, 0.01), { 1.7e-8, 1.8e-6, 39 }));
    CHECK(accurate<4>(pole, 0.0, pole_n(4, 0.01), { 2.3e-8, 3.3e-6, 109 }));

    // float: sin rounds alike in numpy and the C++ libraries, exp does not, so exp is checked for coverage only.
    auto const sinf = [](float t) { return std::sin(t); };
    auto const expf = [](float t) { return std::exp(t); };
    CHECK(accurate<2>(sinf, 0.5f, float(sin_n(2, 0.5)), { 2.5e-5, 1.9e-4, 11 }));
    CHECK(accurate<3>(sinf, 0.5f, float(sin_n(3, 0.5)), { 1.8e-4, 4.7e-3, 20 }));
    CHECK(covered<2>(expf, 1.0f, float(e), 6.0f));
    CHECK(covered<3>(expf, 1.0f, float(e), 6.0f));
    // f'''' in float is the low-precision case: deliberate here, so the warning is silenced at the call.
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(push)
#pragma warning(disable : 4996)
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    auto const noisy = derivative<Ridders, 4>(sinf, 0.5f);
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif
    // The prototype: 1.6e-3 relative off, with an estimate of 5.2e-2 relative; the estimate shows the noise (at least
    // half the prototype's).
    REQUIRE(noisy.has_value());
    CHECK(std::abs(noisy->value - float(sin_n(4, 0.5))) <= 6.0f * *noisy->abs_error);
    CHECK(*noisy->abs_error >= 2.6e-2f * float(sin_n(4, 0.5)));
}

TEST_CASE("orders: stencils with the default step")
{
    // For exp every derivative equals f, so with h = eps^(1/(p+N)) the relative error is about C eps^(p/(p+N)) from
    // truncation plus S eps^(p/(p+N)) from rounding fn's values and up to S/2 more from rounding the points (the
    // nodes stay within |x| + K h <= 1.25 here): (C + 2S) eps^(p/(p+N)), with C = |M_(N+p)| / ((N+p)! den) and S the
    // sum of the absolute weights over den, both computed from the moments here.
    double const eps = std::numeric_limits<double>::epsilon();
    auto const check = [eps]<typename S, int N>(S, std::integral_constant<int, N>) {
        constexpr auto t = table<S, N>;
        constexpr int p = t.order;
        double const c = std::abs(double(moment(t, N + p))) / (double(factorial(N + p)) * t.denominator);
        double s = 0;
        for (int const w : t.numerators) s += std::abs(w);
        s /= t.denominator;
        double const tol = (c + 2 * s) * std::pow(eps, double(p) / (p + N));
        auto const exp = [](double x) { return std::exp(x); };
        for (double const x : { -1.0, 0.5, 1.0 }) {
            auto const d = derivative<S, N>(exp, x);
            REQUIRE(d.has_value());
            CHECK(std::abs(d->value - std::exp(x)) <= tol * std::exp(x));
        }
    };
    using two = std::integral_constant<int, 2>;
    using three = std::integral_constant<int, 3>;
    using four = std::integral_constant<int, 4>;
    check(Central<2>{}, two{});
    check(Central<4>{}, two{});
    check(Central<10>{}, two{});
    check(Forward<1>{}, two{});
    check(Forward<6>{}, two{});
    check(Backward<3>{}, two{});
    check(Central<2>{}, three{});
    check(Central<6>{}, three{});
    check(Forward<2>{}, three{});
    check(Backward<6>{}, three{});
    check(Central<2>{}, four{});
    check(Central<8>{}, four{});
    check(Forward<4>{}, four{});
    check(Backward<1>{}, four{});
}

TEST_CASE("orders: multiprecision types")
{
    // Ridders f'' to f'''' of exp at 1: coverage as for double (the actual error within 6 times the estimate), and
    // the estimate bounded per type by twice the first C++ run, as for f' (docs/design.md: the decimal prototype does
    // not model these types' guard digits). The first run gave relative estimates of 2.7e-43, 1.6e-39 and 8.1e-36 for
    // cpp_bin_float_50 (with and without expression templates alike), 1.0e-42, 3.9e-39 and 2.6e-36 for
    // cpp_dec_float_50, and relative errors 13 to 2600 times smaller.
    auto const check = []<typename T>(T, double bound2, double bound3, double bound4) {
        auto const exp = [](T const &x) { return T(mp::exp(x)); };
        T const e = T(mp::exp(T(1)));
        auto const order_case = [&](auto n, double bound) {
            auto const r = derivative<Ridders, decltype(n)::value>(exp, T(1));
            REQUIRE(r.has_value());
            T const actual = T(mp::abs(T(r->value - e)));
            CHECK(actual <= T(6) * *r->abs_error);
            CHECK(*r->abs_error <= T(bound) * e);
        };
        order_case(std::integral_constant<int, 2>{}, bound2);
        order_case(std::integral_constant<int, 3>{}, bound3);
        order_case(std::integral_constant<int, 4>{}, bound4);
    };
    check(bin50{}, 6e-43, 4e-39, 2e-35);
    check(bin50_et{}, 6e-43, 4e-39, 2e-35);
    check(dec50{}, 3e-42, 8e-39, 6e-36);

    // The floor saturates in these types too, with inputs from the type: a constant max / 4 with h = eps at 0. The
    // first division gives (2 eps max / 4) / eps = max / 2 on each side, their sum max, and the last division by
    // h = eps then saturates, so the result is {0, max}.
    auto const saturates = []<typename T>(T) {
        T const largest = (std::numeric_limits<T>::max)();
        T const eps = std::numeric_limits<T>::epsilon();
        auto const flat = [largest](T const &) { return T(largest / T(4)); };
        auto const r = derivative<Ridders, 2>(flat, T(0), eps);
        REQUIRE(r.has_value());
        CHECK(r->value == T(0));
        CHECK(*r->abs_error == largest);
    };
    saturates(bin50{});
    saturates(dec50{});
}
