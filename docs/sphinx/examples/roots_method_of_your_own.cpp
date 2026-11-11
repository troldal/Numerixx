// The program of docs/sphinx/roots/method-of-your-own.md: the Pegasus method (Dowell and Jarratt 1972) on a bracket
// and Chebyshev's method from a guess, written as stepped methods against the public protocol only, checked against
// their concepts where they are defined, and called through solve() on Wallis's cubic x^3 - 2x - 5. The page names the
// Pegasus method without its source, which the manual's citation rule does not admit (PHASE10_PLAN.md, Q6). The page
// shows the parts between the markers and the output in the raw strings; main returns non-zero unless the program
// prints that output and what the page says holds. The methods and the function use only + - * /, which IEEE 754
// rounds correctly, so the output is the same on every target that does not contract a * b + c into a fused
// multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <array>
#include <cmath>
#include <expected>
#include <functional>
#include <optional>
#include <print>
#include <tuple>

namespace rt = nxx::roots;
// </includes>
// clang-format on

#include "printed.hpp"

#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <string_view>

// clang-format off
// <pegasus>
// The Pegasus method: false position whose retained end's value is scaled by f1 / (f1 + fx) whenever that end is
// retained again; each new point kept `inside` times tol(point) away from the ends
struct Pegasus
{
    double inside = 1.0;

    template<nxx::real T>
    struct State
    {
        T x0, f0, g0; // the retained end, fn's value there, and the scaled value the secant uses
        T x1, f1;     // the latest point and fn's value there
        rt::Step<T> step;
    };

    // The report: both points in order with fn's true values there, the point the one with the smaller |f|
    template<nxx::real T>
    static constexpr rt::Step<T> report(State<T> const &s)
    {
        T const point = nxx::math::abs(s.f1) <= nxx::math::abs(s.f0) ? s.x1 : s.x0;
        if (s.x0 < s.x1) return rt::Step<T>{ s.x0, s.f0, s.x1, s.f1, point };
        return rt::Step<T>{ s.x1, s.f1, s.x0, s.f0, point };
    }

    template<nxx::real T, typename F>
    constexpr std::expected<State<T>, rt::error>
        start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
    {
        State<T> state{ s.lower, s.at_lower[0], s.at_lower[0], s.upper, s.at_upper[0], {} };
        state.step = report(state);
        return state;
    }

    template<nxx::real T, typename F>
    constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &tol) const
    {
        T const lo = s.step.lower;
        T const hi = s.step.upper;
        // The secant through (x1, f1) and (x0, g0), from halves, so that no difference overflows: f1 and g0
        // have opposite signs, so q lies in [0, 1]; d is half the step, and x = x1 - 2d lies between the ends,
        // formed in two subtractions where 2d itself would overflow (a bracket wider than max)
        T const half1 = T(s.f1 / T(2));
        T const q = T(half1 / T(half1 - T(s.g0 / T(2))));
        T const d = T(q * T(T(s.x1 / T(2)) - T(s.x0 / T(2))));
        T const quarter = T(T((std::numeric_limits<T>::max)()) / T(4));
        T x = nxx::math::abs(d) <= quarter ? T(s.x1 - T(d + d)) : T(T(s.x1 - d) - d);
        T const margin = T(T(inside) * tol.at(s.step.point));
        if (x < T(lo + margin)) x = T(lo + margin);
        else if (T(hi - margin) < x) x = T(hi - margin);
        // The midpoint where the bracket is too narrow for the margin or x is not inside; halves again
        if (T(T(hi / T(4)) - T(lo / T(4))) <= margin || !(lo < x && x < hi)) x = T(T(lo / T(2)) + T(hi / T(2)));
        if (!(lo < x && x < hi)) return std::unexpected(rt::error::not_converged);
        T const fx = T(std::invoke(fn, x));
        if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
        if (fx == T(0)) {
            s.step = rt::Step<T>{ lo, s.step.f_lower, x, fx, x };
            return {};
        }
        if ((fx < T(0)) != (s.f1 < T(0))) { // the sign test without a product, which could underflow
            s.x0 = s.x1;
            s.f0 = s.f1;
            s.g0 = s.f1;
        } else { // f1 and fx have one sign: halved, their sum cannot overflow
            s.g0 = T(s.g0 * T(half1 / T(half1 + T(fx / T(2)))));
        }
        s.x1 = x;
        s.f1 = fx;
        s.step = report(s);
        return {};
    }
};
static_assert(rt::bracketing_method<Pegasus, double>); // checked where it is defined
// </pegasus>
// clang-format on

// "Nothing overflows": in constant evaluation an overflow is not a constant, so these compile only if none occurs.
// Values near max (the independent review's case: 1.6e308 (2x - 3) on [1, 2], root 1.5, where f1 - g0 overflowed
// before the halves), and a bracket wider than max (x - 1.5 on [-1e308, 1e308], where x1 - x0 would).
namespace {

    constexpr bool brackets_root(std::expected<rt::Root<double>, rt::error> const &r, double root)
    { return r.has_value() && nxx::math::abs(r->value - root) <= r->abs_error; }

    static_assert(brackets_root(
        rt::solve<Pegasus>([](double x) { return 1.6e308 * (2.0 * x - 3.0); }, rt::Bracket{ 1.0, 2.0 }), 1.5));
    static_assert(
        brackets_root(rt::solve<Pegasus>([](double x) { return x - 1.5; }, rt::Bracket{ -1e308, 1e308 }), 1.5));

} // namespace

// clang-format off
// <chebyshev>
// Chebyshev's method (values = 3): the step -t (1 + t f'' / (2 f')), t = f / f', at least `min_step` times
// tol(x) long; at the first sign change between two iterates it continues as Brent's method, whose state it holds
struct Chebyshev
{
    static constexpr int values = 3;
    double min_step = 1.0;

    template<nxx::real T>
    struct State
    {
        T x;                // the latest iterate
        std::array<T, 3> v; // f, f' and f'' there
        std::optional<rt::bracketing_state_t<rt::Brent, T>> brent; // after the first sign change
        rt::Step<T> step;
    };

    template<nxx::real T, typename F>
    static constexpr std::array<T, 3> values_at(F &fn, T const &x)
    {
        auto const r = std::invoke(fn, x);
        return { T(std::get<0>(r)), T(std::get<1>(r)), T(std::get<2>(r)) };
    }

    template<nxx::real T, typename F>
    constexpr std::expected<State<T>, rt::error> start(F &fn, T const &guess, rt::Tolerance<T> const &) const
    {
        std::array<T, 3> const v = values_at(fn, guess);
        if (!nxx::math::isfinite(v[0])) return std::unexpected(rt::error::non_finite);
        return State<T>{ guess, v, std::nullopt, rt::Step<T>{ guess, v[0], guess, v[0], guess } };
    }

    template<nxx::real T, typename F>
    constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &tol) const
    {
        if (s.brent) { // Brent's step, and its report
            std::expected<void, rt::error> const r = rt::Brent{}.advance(fn, *s.brent, tol);
            s.step = s.brent->step;
            return r;
        }
        if (s.v[1] == T(0) || !nxx::math::isfinite(s.v[1]) || !nxx::math::isfinite(s.v[2]))
            return std::unexpected(rt::error::not_converged);
        T const t = T(s.v[0] / s.v[1]);
        T step = T(T(-t) * T(T(1) + T(T(t * s.v[2]) / T(T(2) * s.v[1]))));
        T const least = T(T(min_step) * tol.at(s.x));
        if (nxx::math::abs(step) < least) step = step < T(0) ? T(-least) : least;
        T const next = T(s.x + step);
        std::array<T, 3> const vn = values_at(fn, next);
        if (!nxx::math::isfinite(vn[0])) return std::unexpected(rt::error::non_finite);
        bool const up = s.x < next;
        if (vn[0] != T(0) && (vn[0] < T(0)) != (s.v[0] < T(0))) { // a sign change: Brent from here
            rt::Start<T> const ends = up ? rt::Start<T>{ s.x, { s.v[0] }, next, { vn[0] } }
                                         : rt::Start<T>{ next, { vn[0] }, s.x, { s.v[0] } };
            auto started = rt::Brent{}.start(fn, ends, tol);
            if (!started) return std::unexpected(started.error());
            s.brent = *started;
            s.step = s.brent->step;
            return {};
        }
        s.step = up ? rt::Step<T>{ s.x, s.v[0], next, vn[0], next }
                    : rt::Step<T>{ next, vn[0], s.x, s.v[0], next };
        s.x = next;
        s.v = vn;
        return {};
    }
};
static_assert(rt::open_method<Chebyshev, double>);
// </chebyshev>
// clang-format on

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <call>
    int n = 0; // evaluations of the cubic
    auto const cubic = [&n](double x) { ++n; return (x * x - 2.0) * x - 5.0; };
    auto const print = [&n](char const *name, std::expected<rt::Root<double>, rt::error> const &r) {
        if (r) std::println("{:16} {:.17g} +- {:.1e}, {} evaluations", name, r->value, r->abs_error, n);
        else std::println("{:16} {}, {} evaluations", name, r.error(), n);
        n = 0;
    };
    print("Pegasus", rt::solve<Pegasus>(cubic, rt::Bracket{ 2.0, 3.0 }));
    print("inside = 0", rt::solve(Pegasus{ .inside = 0.0 }, cubic, rt::Bracket{ 2.0, 3.0 }));
    print("AndersonBjorck", rt::solve<rt::AndersonBjorck>(cubic, rt::Bracket{ 2.0, 3.0 }));
    // </call>
    // clang-format on

    constexpr std::string_view call_output = R"(
Pegasus          2.0945514815423265 +- 8.9e-16, 10 evaluations
inside = 0       2.0945514815423265 +- 1.3e-15, 33 evaluations
AndersonBjorck   2.0945514815423265 +- 1.3e-15, 9 evaluations
)";

    // clang-format off
    // <open-call>
    auto const df = [](double x) { return 3.0 * x * x - 2.0; };
    auto const d2f = [](double x) { return 6.0 * x; };
    print("Chebyshev", rt::solve<Chebyshev>(cubic, df, d2f, 2.0));
    auto const all = [&cubic](double x) { return std::array{ cubic(x), 3.0 * x * x - 2.0, 6.0 * x }; };
    print("min_step = 0", rt::solve(Chebyshev{ .min_step = 0.0 }, all, 2.0));
    print("Halley", rt::solve<rt::Halley>(cubic, df, d2f, 2.0));
    // </open-call>
    // clang-format on

    constexpr std::string_view open_output = R"(
Chebyshev        2.0945514815423265 +- 8.9e-16, 5 evaluations
min_step = 0     not_converged, 213 evaluations
Halley           2.0945514815423265 +- 8.9e-16, 5 evaluations
)";

    // clang-format off
    // <observe>
    auto const show = [](rt::Step<double> const &s) {
        std::println("[{:.17g}, {:.17g}] point {:.17g}", s.lower, s.upper, s.point);
    };
    auto const r = rt::solve<Pegasus>(rt::observer(show), cubic, rt::Bracket{ 2.0, 3.0 }, { .rel_tol = 1e-6 });
    // </observe>
    // clang-format on
    n = 0;

    constexpr std::string_view observe_output = R"(
[2, 3] point 2
[2.0588235294117645, 3] point 2.0588235294117645
[2.0812636598450229, 3] point 2.0812636598450229
[2.0927546010136555, 3] point 2.0927546010136555
[2.0945279947670867, 3] point 2.0945279947670867
[2.0945279947670867, 2.0945515006841915] point 2.0945515006841915
[2.0945494061326908, 2.0945515006841915] point 2.0945515006841915
)";

    std::string expected;
    expected += docs::quoted(call_output);
    expected += docs::quoted(open_output);
    expected += docs::quoted(observe_output);
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    // "Each result lies within abs_error of Wallis's root": 2.0945514815423265915... (OEIS A007493); the cubic's
    // rounding near the root, about 10 eps of its terms of up to 9 over |f'| = 11.2, moves its sign change by under
    // 1e-15, below abs_error plus that.
    double const wallis = 2.0945514815423265915;
    double const eps = std::numeric_limits<double>::epsilon();
    auto const near = [&](std::expected<rt::Root<double>, rt::error> const &x) {
        return x && std::abs(x->value - wallis) <= x->abs_error + 10.0 * eps * 9.0 / 11.2;
    };
    auto const plain = [](double x) { return (x * x - 2.0) * x - 5.0; };
    if (!(near(rt::solve<Pegasus>(plain, rt::Bracket{ 2.0, 3.0 }))
          && near(rt::solve(Pegasus{ .inside = 0.0 }, plain, rt::Bracket{ 2.0, 3.0 }))
          && near(rt::solve<Chebyshev>(plain, df, d2f, 2.0)) && near(r)))
        return 2;
    // "Without its minimum step, from 2 its iterates approach the root from below and never cross it": fn is
    // negative at every point reported, and the iterates never decrease (they end on the double below the root,
    // where the step no longer moves them).
    bool below = true;
    double last = 0.0;
    auto const watch = [&](rt::Step<double> const &s) {
        below = below && s.f_lower < 0.0 && s.f_upper < 0.0 && last <= s.point;
        last = s.point;
    };
    if (!(rt::solve(Chebyshev{ .min_step = 0.0 }, rt::observer(watch), all, 2.0)
              == std::unexpected(rt::error::not_converged)
          && below))
        return 3;
    // "The upper end stays at 3 for the start and the four steps after it": the observed run's first five reports.
    int at_three = 0;
    auto const count = [&at_three](rt::Step<double> const &s) { at_three += s.upper == 3.0 ? 1 : 0; };
    (void)rt::solve<Pegasus>(rt::observer(count), plain, rt::Bracket{ 2.0, 3.0 }, { .rel_tol = 1e-6 });
    if (!(at_three == 5)) return 4;
    // "An observed solve computes exactly what an unobserved one does".
    if (!(r == rt::solve<Pegasus>(plain, rt::Bracket{ 2.0, 3.0 }, { .rel_tol = 1e-6 }))) return 5;
    // The reports of a run, for the claims about the stall and the crossing.
    auto const reports = [&plain](Pegasus const &m, rt::Criteria<double> const &c) {
        std::vector<rt::Step<double>> steps;
        auto record = [&steps](rt::Step<double> const &s) { steps.push_back(s); };
        (void)rt::solve(m, rt::observer(record), plain, rt::Bracket{ 2.0, 3.0 }, c);
        return steps;
    };
    // "With inside = 0, after nine evaluations its lower end is the double just below the root, and only the
    // midpoint fallback narrows the bracket, by half, for the remaining 24": the report after 7 steps (2 + 7
    // evaluations) has that lower end, fn negative there and positive at the next double; it never moves again, and
    // each of the 24 later reports is at most half as wide as the one before, plus one unit at the root (the last
    // halvings round: 5 units to 3).
    std::vector<rt::Step<double>> const stall = reports(Pegasus{ .inside = 0.0 }, {});
    if (!(stall.size() == 32)) return 6;
    double const stuck = stall[7].lower;
    if (!(plain(stuck) < 0.0 && plain(std::nextafter(stuck, 3.0)) > 0.0 && stall[6].lower != stuck)) return 7;
    for (std::size_t i = 8; i < stall.size(); ++i) {
        double const width = stall[i].upper - stall[i].lower;
        double const before = stall[i - 1].upper - stall[i - 1].lower;
        if (!(stall[i].lower == stuck && width <= 0.5 * before + 2.0 * eps * stuck)) return 8;
    }
    // "The scaled value tilts the line, and the point lands above the root; the last point, kept tol inside the
    // upper end, closes the bracket": the observed run's first six reports are those of inside = 0, the sixth's point
    // above the root (fn positive), and its last lower end is the upper end less tol there.
    std::vector<rt::Step<double>> const seen = reports(Pegasus{}, { .rel_tol = 1e-6 });
    std::vector<rt::Step<double>> const plain_seen = reports(Pegasus{ .inside = 0.0 }, { .rel_tol = 1e-6 });
    if (!(seen.size() == 7 && plain_seen.size() >= 6)) return 9;
    for (std::size_t i = 0; i < 6; ++i)
        if (!(seen[i] == plain_seen[i])) return 10;
    rt::Tolerance<double> const loose{ std::numeric_limits<double>::denorm_min(), 1e-6 };
    if (!(plain(seen[5].point) > 0.0 && seen[6].lower == seen[5].upper - loose.at(seen[5].point))) return 11;
    // The concepts: each method is of its kind only.
    static_assert(!rt::open_method<Pegasus, double> && !rt::bracketing_method<Chebyshev, double>);
    static_assert(rt::values_of<Chebyshev> == 3 && rt::values_of<Pegasus> == 1);
    return 0;
}
