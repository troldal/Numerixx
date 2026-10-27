// The one-sided extrapolation methods, ForwardRidders and BackwardRidders (Phase 5, docs/deriv/PHASE5_API.md).
//
// References are closed forms, and the independent Python prototype .scratch/onesided/ridders5.py (numpy, written
// from the textbook rules), which reproduces the C++ results for double exactly on 28,000 cases: values, estimates
// and evaluations. Bounds as in test_ridders.cpp and test_orders.cpp: the actual error within 6 times the estimate
// (the prototype's worst over 10,800 random cases in double was 2.6 times), within 4 times the prototype's error, the
// estimate within twice the prototype's, and the evaluations within the prototype's plus 4 (libm implementations
// round fn's values differently, which can move a stall by a level).

#include <numerixx/deriv.hpp>

#include "silence.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <expected>
#include <limits>
#include <numbers>
#include <type_traits>

namespace mp = boost::multiprecision;

using bin50 = mp::cpp_bin_float_50;
using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
using dec50 = mp::cpp_dec_float_50;

using nxx::deriv::BackwardRidders;
using nxx::deriv::derivative;
using nxx::deriv::Estimate;
using nxx::deriv::ForwardRidders;
using nxx::deriv::order;
using nxx::deriv::Ridders;

namespace {

    template<int K>
    constexpr auto power = [](auto x) {
        decltype(x) p = 1;
        for (int i = 0; i < K; ++i) p = p * x;
        return p;
    };

    constexpr long long factorial(int n) { return n <= 1 ? 1 : n * factorial(n - 1); }

    constexpr auto non_finite = std::unexpected(nxx::deriv::error::non_finite);

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

    // Records the smallest and the largest point fn is evaluated at.
    template<typename F>
    struct Recorded
    {
        F fn;
        double *lowest;
        double *highest;
        double operator()(double x) const
        {
            *lowest = std::min(*lowest, x);
            *highest = std::max(*highest, x);
            return fn(x);
        }
    };

    // The prototype's values for a case, relative to |f^(N)|, and its number of evaluations.
    struct Reference
    {
        double error;
        double estimate;
        int evaluations;
    };

    template<typename T>
    bool judged(std::expected<Estimate<T>, nxx::deriv::error> const &r, int calls, T exact, Reference ref)
    {
        if (!r || !r->abs_error) return false;
        T const actual = std::abs(r->value - exact);
        bool const ok = actual <= T(6) * *r->abs_error && actual <= T(4) * T(ref.error) * std::abs(exact)
                        && *r->abs_error <= T(2) * T(ref.estimate) * std::abs(exact) && calls <= ref.evaluations + 4;
        if (!ok) MESSAGE("actual ", actual, ", estimate ", *r->abs_error, ", calls ", calls);
        return ok;
    }

    template<int N, typename F>
    bool accurate(F fn, double x, double exact, Reference ref)
    {
        int calls = 0;
        auto const r = derivative<ForwardRidders, N>(Counted{ fn, &calls }, x);
        return judged(r, calls, exact, ref);
    }

    template<typename T>
    bool covered(std::expected<Estimate<T>, nxx::deriv::error> const &r, T exact)
    { return r && r->abs_error && std::abs(r->value - exact) <= T(6) * *r->abs_error; }

    // Closed forms of the N-th derivative.
    double sin_n(int n, double x) { return std::sin(x + n * std::numbers::pi / 2); }
    double log_n(int n, double x) { return (n % 2 == 1 ? 1 : -1) * double(factorial(n - 1)) / std::pow(x, n); }
    double sqrt_n(int n, double x)
    {
        double c = 1;
        for (int k = 0; k < n; ++k) c *= 0.5 - k;
        return c * std::pow(x, 0.5 - n);
    }
    double pole_n(int n, double a) { return -double(factorial(n)) / std::pow(a, n + 1); } // 1/(x - a) at 0
    // 1/(1 + 25 x^2) = Re 1/(1 + 5ix), whose N-th derivative is (-1)^N N! (5i)^N (1 + 5ix)^-(N+1).
    double runge_n(int n, double x)
    {
        std::complex<double> const z(1.0, 5.0 * x);
        std::complex<double> const c = std::pow(std::complex<double>(0.0, 5.0), n) * double(factorial(n));
        return ((n % 2 == 0 ? c : -c) / std::pow(z, n + 1)).real();
    }

    // The forward rule of order N is exact for x^N at every step, so the first two levels agree exactly: N! at x = 2
    // and h = 1/8, after 1 + 2N evaluations (fn(x), then N per level). x^(N+1) has the error (N+1)! N/2 h, linear in
    // h, which the first column removes exactly: 2 (N+1)! after three levels, 1 + 3N evaluations. Every value is a
    // short dyadic fraction, exact in double, and so are the divisions by powers of h and by 2^j - 1 here.
    template<int N, typename M>
    constexpr bool exact_on_polynomials()
    {
        auto const low = derivative<M, N>(power<N>, 2.0, 0.125);
        auto const high = derivative<M, N>(power<N + 1>, 2.0, 0.125);
        return low && low->value == double(factorial(N)) && high && high->value == 2.0 * double(factorial(N + 1));
    }

    // The estimate of an exact stop is the floor carried through the tableau: for f' of x at 2 the levels' floors are
    // F_i = (eps |f(x)| + eps |f(x + h_i)|) / h_i, and the entry A[1][1] = 2 D_1 - D_0 has the floor 2 F_1 + F_0.
    constexpr bool floor_of_exact_stop()
    {
        double const eps = std::numeric_limits<double>::epsilon();
        auto const f0 = [eps](double h) { return (eps * 2.0 + eps * (2.0 + h)) / h; };
        auto const r = derivative<ForwardRidders>(power<1>, 2.0, 0.125);
        return r && r->value == 1.0 && *r->abs_error == 2.0 * f0(0.0625) + f0(0.125);
    }

    // A kink at 1: t^2 below, 3t - 2 above. The one-sided methods return the one-sided derivatives, exactly: above, the
    // rule is exact for a line; below, the error -h of the forward rule on t^2 is removed by the first column.
    constexpr auto kink = [](double t) { return t < 1.0 ? t * t : 3.0 * t - 2.0; };

} // namespace

// The new methods are accepted for the same orders as Ridders.
template<typename M, auto N>
concept accepted = requires { derivative<M, N>(power<2>, 1.0); };
static_assert(accepted<ForwardRidders, 1> && accepted<ForwardRidders, order::fourth> && accepted<BackwardRidders, 3>);
static_assert(!accepted<ForwardRidders, 0> && !accepted<BackwardRidders, 5>);

static_assert(exact_on_polynomials<1, ForwardRidders>() && exact_on_polynomials<2, ForwardRidders>());
static_assert(exact_on_polynomials<3, ForwardRidders>() && exact_on_polynomials<4, ForwardRidders>());
static_assert(exact_on_polynomials<1, BackwardRidders>() && exact_on_polynomials<2, BackwardRidders>());
static_assert(exact_on_polynomials<3, BackwardRidders>() && exact_on_polynomials<4, BackwardRidders>());
static_assert(floor_of_exact_stop());
static_assert(derivative<ForwardRidders>(kink, 1.0, 1.0 / 32)->value == 3.0);
static_assert(derivative<BackwardRidders>(kink, 1.0, 1.0 / 32)->value == 2.0);
// float, where the values stay exact: x^4 at 1 with h = 1/8 and 1/16 (at most (20/16)^4 = 160000 / 2^16, below 2^24).
// Deliberate: a float f'''' warns, so it is silenced.
NUMERIXX_SILENCE_BEGIN
static_assert(derivative<ForwardRidders, 4>(power<4>, 1.0f, 0.125f)->value == 24.0f);
NUMERIXX_SILENCE_END

TEST_CASE("onesided: exact on polynomials, with the evaluations")
{
    CHECK(exact_on_polynomials<1, ForwardRidders>());
    CHECK(exact_on_polynomials<4, ForwardRidders>());
    CHECK(exact_on_polynomials<2, BackwardRidders>());
    CHECK(exact_on_polynomials<3, BackwardRidders>());
    CHECK(floor_of_exact_stop());
    auto const evaluations = []<int N>(std::integral_constant<int, N>) {
        int low = 0;
        int high = 0;
        auto const a = derivative<ForwardRidders, N>(Counted{ power<N>, &low }, 2.0, 0.125);
        auto const b = derivative<ForwardRidders, N>(Counted{ power<N + 1>, &high }, 2.0, 0.125);
        return a && b && low == 1 + 2 * N && high == 1 + 3 * N;
    };
    CHECK(evaluations(std::integral_constant<int, 1>{}));
    CHECK(evaluations(std::integral_constant<int, 2>{}));
    CHECK(evaluations(std::integral_constant<int, 3>{}));
    CHECK(evaluations(std::integral_constant<int, 4>{}));
    CHECK(derivative<ForwardRidders>(kink, 1.0)->value == 3.0);
    CHECK(derivative<BackwardRidders>(kink, 1.0)->value == 2.0);
}

TEST_CASE("onesided: accuracy against closed forms, double and float")
{
    auto const exp = [](double t) { return std::exp(t); };
    auto const sin = [](double t) { return std::sin(t); };
    auto const log = [](double t) { return std::log(t); };
    auto const runge = [](double t) { return 1.0 / (1.0 + 25.0 * t * t); };
    double const e = std::numbers::e;

    // Prototype references {relative error, relative estimate, evaluations}, from the default initial step s/32.
    CHECK(accurate<1>(exp, 1.0, e, { 7.1e-14, 2.4e-12, 9 }));
    CHECK(accurate<1>(sin, 0.5, sin_n(1, 0.5), { 3.3e-14, 1.3e-12, 9 }));
    CHECK(accurate<1>(runge, 0.3, runge_n(1, 0.3), { 5.7e-14, 2.1e-12, 11 }));
    CHECK(accurate<1>(log, 0.2, log_n(1, 0.2), { 2.3e-15, 3.1e-12, 11 }));
    CHECK(accurate<2>(exp, 1.0, e, { 1.1e-11, 1.3e-9, 15 }));
    CHECK(accurate<2>(sin, 0.5, sin_n(2, 0.5), { 2.4e-11, 1.3e-9, 15 }));
    CHECK(accurate<2>(runge, 0.3, runge_n(2, 0.3), { 1.3e-10, 1.2e-9, 19 }));
    CHECK(accurate<2>(log, 0.2, log_n(2, 0.2), { 3.6e-11, 4.1e-9, 21 }));
    CHECK(accurate<3>(exp, 1.0, e, { 1.4e-7, 1.4e-6, 22 }));
    CHECK(accurate<3>(sin, 0.5, sin_n(3, 0.5), { 7.6e-8, 1.8e-6, 22 }));
    CHECK(accurate<3>(runge, 0.3, runge_n(3, 0.3), { 3.6e-8, 1.4e-6, 28 }));
    CHECK(accurate<3>(log, 0.2, log_n(3, 0.2), { 1.7e-8, 1.1e-6, 28 }));
    CHECK(accurate<4>(exp, 1.0, e, { 4.7e-7, 2.2e-3, 25 }));
    CHECK(accurate<4>(sin, 0.5, sin_n(4, 0.5), { 4.9e-6, 4.9e-4, 25 }));
    CHECK(accurate<4>(runge, 0.3, runge_n(4, 0.3), { 5.8e-5, 3.8e-3, 37 }));
    CHECK(accurate<4>(log, 0.2, log_n(4, 0.2), { 3.5e-6, 1.1e-4, 33 }));

    // float: accuracy for f' (sin rounds alike in numpy and the C++ libraries, exp does not, so exp is checked for
    // coverage only), coverage for f'' (worst 4.0 times the estimate in the prototype's random study). The higher
    // orders in float are the noise the low-precision warning is about; f'' warns too, so it is silenced here.
    auto const sinf = [](float t) { return std::sin(t); };
    auto const expf = [](float t) { return std::exp(t); };
    int calls = 0;
    auto const sin1 = derivative<ForwardRidders>(Counted{ sinf, &calls }, 0.5f);
    CHECK(judged(sin1, calls, float(sin_n(1, 0.5)), { 5.8e-6, 8.4e-5, 6 }));
    CHECK(covered(derivative<ForwardRidders>(expf, 1.0f), float(e)));
    NUMERIXX_SILENCE_BEGIN
    auto const sin2 = derivative<ForwardRidders, 2>(sinf, 0.5f);
    auto const exp2 = derivative<ForwardRidders, 2>(expf, 1.0f);
    NUMERIXX_SILENCE_END
    CHECK(covered(sin2, float(sin_n(2, 0.5))));
    CHECK(covered(exp2, float(e)));
}

TEST_CASE("onesided: domain edges")
{
    // fn undefined below x = 1: every central sample reaches below it, the one-sided ones never do.
    auto const exp_above = [](double t) { return t < 1.0 ? std::numeric_limits<double>::quiet_NaN() : std::exp(t); };
    auto const exp_below = [](double t) { return t > 1.0 ? std::numeric_limits<double>::quiet_NaN() : std::exp(t); };
    double const e = std::numbers::e;
    CHECK(derivative<Ridders>(exp_above, 1.0) == non_finite);
    CHECK(accurate<1>(exp_above, 1.0, e, { 7.1e-14, 2.4e-12, 9 }));
    CHECK(accurate<2>(exp_above, 1.0, e, { 1.1e-11, 1.3e-9, 15 }));
    CHECK(accurate<3>(exp_above, 1.0, e, { 1.4e-7, 1.4e-6, 22 }));
    CHECK(accurate<4>(exp_above, 1.0, e, { 4.7e-7, 2.2e-3, 25 }));
    CHECK(covered(derivative<BackwardRidders>(exp_below, 1.0), e));
    CHECK(covered(derivative<BackwardRidders, 4>(exp_below, 1.0), e));
    // The example of ForwardRidders' comment in deriv.hpp: log, NaN below 0.5, at 0.5. Ridders fails, ForwardRidders
    // gives 2 to about 1e-13 relative (4 * 9.4e-15 here; the prototype's run of this case gives the reference).
    auto const log_above = [](double t) { return t < 0.5 ? std::nan("") : std::log(t); };
    CHECK(derivative<Ridders>(log_above, 0.5) == non_finite);
    CHECK(accurate<1>(log_above, 0.5, 2.0, { 9.4e-15, 1.7e-12, 10 }));

    // sqrt and log near 0, where fn varies on the scale x, far below the default step: after untrusted results the
    // restarts start below the smallest step tried, and reach the scale of x within the levels.
    auto const sqrt = [](double t) { return std::sqrt(t); };
    auto const log = [](double t) { return std::log(t); };
    CHECK(accurate<1>(sqrt, 1e-4, sqrt_n(1, 1e-4), { 1.3e-13, 2.0e-12, 21 }));
    CHECK(accurate<1>(log, 1e-4, log_n(1, 1e-4), { 6.0e-13, 9.4e-12, 21 }));
    CHECK(accurate<1>(sqrt, 1e-8, sqrt_n(1, 1e-8), { 6.4e-14, 1.7e-12, 33 }));
    CHECK(accurate<1>(log, 1e-8, log_n(1, 1e-8), { 1.8e-12, 1.5e-11, 33 }));
    CHECK(accurate<2>(sqrt, 1e-4, sqrt_n(2, 1e-4), { 1.8e-10, 8.7e-9, 41 }));
    CHECK(accurate<2>(log, 1e-4, log_n(2, 1e-4), { 1.7e-11, 1.4e-8, 41 }));
    CHECK(accurate<2>(sqrt, 1e-8, sqrt_n(2, 1e-8), { 6.7e-11, 1.8e-9, 65 }));
    CHECK(accurate<2>(log, 1e-8, log_n(2, 1e-8), { 5.8e-10, 2.3e-8, 65 }));
    CHECK(accurate<3>(sqrt, 1e-4, sqrt_n(3, 1e-4), { 8.6e-8, 1.5e-5, 61 }));
    CHECK(accurate<3>(log, 1e-4, log_n(3, 1e-4), { 6.3e-8, 2.3e-5, 61 }));
    CHECK(accurate<3>(sqrt, 1e-8, sqrt_n(3, 1e-8), { 6.5e-8, 9.1e-7, 97 }));
    CHECK(accurate<3>(log, 1e-8, log_n(3, 1e-8), { 6.3e-7, 2.6e-6, 97 }));
    CHECK(accurate<4>(sqrt, 1e-4, sqrt_n(4, 1e-4), { 1.6e-6, 7.8e-4, 77 }));
    CHECK(accurate<4>(log, 1e-4, log_n(4, 1e-4), { 1.4e-5, 2.5e-4, 77 }));
    CHECK(accurate<4>(sqrt, 1e-8, sqrt_n(4, 1e-8), { 4.9e-6, 1.7e-4, 129 }));
    CHECK(accurate<4>(log, 1e-8, log_n(4, 1e-8), { 2.9e-6, 2.3e-4, 129 }));

    // A pole 1e-4 behind x, never sampled, which makes fn vary on that scale.
    auto const behind = [](double t) { return 1.0 / (t + 1e-4); };
    CHECK(accurate<1>(behind, 0.0, pole_n(1, -1e-4), { 8.4e-14, 2.0e-12, 22 }));
    CHECK(accurate<2>(behind, 0.0, pole_n(2, -1e-4), { 4.6e-11, 1.3e-9, 43 }));
    CHECK(accurate<3>(behind, 0.0, pole_n(3, -1e-4), { 6.6e-9, 2.2e-6, 64 }));
    CHECK(accurate<4>(behind, 0.0, pole_n(4, -1e-4), { 8.4e-7, 6.8e-4, 81 }));

    // Beyond the reach of the default step, log at 1e-12: the levels run out and the result is wrong (90 to 100%
    // off in the prototype), but flagged by an estimate of at least |value| / 4.
    auto const flagged = []<int N>(std::integral_constant<int, N>, auto fn) {
        int calls = 0;
        auto const r = derivative<ForwardRidders, N>(Counted{ fn, &calls }, 1e-12);
        return r && 4 * *r->abs_error >= std::abs(r->value) && calls == 1 + 32 * N;
    };
    CHECK(flagged(std::integral_constant<int, 1>{}, log));
    CHECK(flagged(std::integral_constant<int, 2>{}, log));
    CHECK(flagged(std::integral_constant<int, 3>{}, log));
    CHECK(flagged(std::integral_constant<int, 4>{}, log));
}

TEST_CASE("onesided: fn is evaluated on one side of x only, restarts included")
{
    auto const sides = []<int N>(std::integral_constant<int, N>, auto fn, double x) {
        double lowest = std::numeric_limits<double>::infinity();
        double highest = -lowest;
        auto const forward = derivative<ForwardRidders, N>(Recorded{ fn, &lowest, &highest }, x);
        bool const above = forward.has_value() && lowest == x && highest > x;
        lowest = std::numeric_limits<double>::infinity();
        highest = -lowest;
        auto const mirrored = [&fn](double t) { return fn(-t); };
        auto const backward = derivative<BackwardRidders, N>(Recorded{ mirrored, &lowest, &highest }, -x);
        return above && backward.has_value() && highest == -x && lowest < -x;
    };
    // log at 1e-6 restarts after untrusted results; the pole 1e-4 behind 0 too.
    auto const log = [](double t) { return std::log(t); };
    auto const behind = [](double t) { return 1.0 / (t + 1e-4); };
    CHECK(sides(std::integral_constant<int, 1>{}, log, 1e-6));
    CHECK(sides(std::integral_constant<int, 2>{}, log, 1e-6));
    CHECK(sides(std::integral_constant<int, 3>{}, log, 1e-6));
    CHECK(sides(std::integral_constant<int, 4>{}, log, 1e-6));
    CHECK(sides(std::integral_constant<int, 1>{}, behind, 0.0));
    CHECK(sides(std::integral_constant<int, 4>{}, behind, 0.0));
}

TEST_CASE("onesided: BackwardRidders is ForwardRidders on the mirrored function")
{
    // The same points, so the same values and estimates, bit for bit; the N-th derivative changes sign for odd N.
    auto const runge = [](double t) { return 1.0 / (1.0 + 25.0 * t * t); };
    auto const mirrored = [&runge](double t) { return runge(-t); };
    auto const same = []<int N>(std::integral_constant<int, N>, auto f, auto g, double x) {
        auto const b = derivative<BackwardRidders, N>(f, x);
        auto const fw = derivative<ForwardRidders, N>(g, -x);
        return b && fw && b->value == (N % 2 == 1 ? -fw->value : fw->value) && b->abs_error == fw->abs_error;
    };
    CHECK(same(std::integral_constant<int, 1>{}, runge, mirrored, 0.3));
    CHECK(same(std::integral_constant<int, 2>{}, runge, mirrored, 0.3));
    CHECK(same(std::integral_constant<int, 3>{}, runge, mirrored, 0.3));
    CHECK(same(std::integral_constant<int, 4>{}, runge, mirrored, 0.3));
    CHECK(same(std::integral_constant<int, 3>{}, runge, mirrored, -2.5));
}

TEST_CASE("onesided: failures and cost")
{
    int calls = 0;
    auto const counted = [&calls](auto f) {
        calls = 0;
        return Counted{ f, &calls };
    };
    auto const identity = [](double t) { return t; };
    double const dmax = std::numeric_limits<double>::max();
    double const qnan = std::numeric_limits<double>::quiet_NaN();

    // fn(x) is part of every level: if it is not finite, no step avoids it, after 1 evaluation.
    auto const hole = [qnan](double t) { return t == 0.0 ? qnan : t; };
    CHECK(derivative<ForwardRidders>(counted(hole), 0.0) == non_finite);
    CHECK(calls == 1);
    CHECK(derivative<BackwardRidders, 4>(counted(hole), 0.0) == non_finite);
    CHECK(calls == 1);

    // fn finite at x only: every level fails, and the cap ends it after 32 levels, with fn(x) once. A level stops at
    // its first non-finite value, except the f'' level, which evaluates both of its points first. (The cap with every
    // evaluation made, 1 + 32 N, is the beyond-reach case above.)
    auto const only_at_1 = [qnan](double t) { return t == 1.0 ? 1.0 : qnan; };
    CHECK(derivative<ForwardRidders, 1>(counted(only_at_1), 1.0) == non_finite);
    CHECK(calls == 33);
    CHECK(derivative<ForwardRidders, 2>(counted(only_at_1), 1.0) == non_finite);
    CHECK(calls == 65);
    CHECK(derivative<BackwardRidders, 3>(counted(only_at_1), 1.0) == non_finite);
    CHECK(calls == 33);
    CHECK(derivative<ForwardRidders, 4>(counted(only_at_1), 1.0) == non_finite);
    CHECK(calls == 33);

    // x = max: every step fails the overflow check, so only fn(x) is evaluated. h below the resolution of x ends the
    // first tableau before its first level. A non-finite x: no evaluation.
    CHECK(derivative<ForwardRidders>(counted(identity), dmax) == non_finite);
    CHECK(calls == 1);
    CHECK(derivative<BackwardRidders>(counted(identity), dmax) == non_finite);
    CHECK(calls == 1);
    CHECK(derivative<ForwardRidders>(counted(identity), 1.0, 1e-300) == non_finite);
    CHECK(calls == 1);
    CHECK(derivative<BackwardRidders>(counted(identity), qnan) == non_finite);
    CHECK(calls == 0);
}

TEST_CASE("onesided: multiprecision types")
{
    // exp at 1, every order: the actual error within 6 times the estimate, and a trusted result (an estimate below
    // |value| / 4, the method's own test). No regression baseline: the prototype does not model these types.
    auto const check = []<typename T>(T) {
        auto const exp = [](T const &x) { return T(mp::exp(x)); };
        T const e = T(mp::exp(T(1)));
        auto const good = [](auto const &r, T const &exact) {
            using std::abs;
            return r && r->abs_error && T(abs(T(r->value - exact))) <= T(6) * *r->abs_error
                   && T(*r->abs_error * T(4)) < T(abs(r->value));
        };
        return good(derivative<ForwardRidders, 1>(exp, T(1)), e) && good(derivative<ForwardRidders, 2>(exp, T(1)), e)
               && good(derivative<ForwardRidders, 3>(exp, T(1)), e) && good(derivative<ForwardRidders, 4>(exp, T(1)), e)
               && good(derivative<BackwardRidders, 1>(exp, T(1)), e)
               && good(derivative<BackwardRidders, 4>(exp, T(1)), e);
    };
    CHECK(check(bin50{}));
    CHECK(check(bin50_et{}));
    CHECK(check(dec50{}));
}
