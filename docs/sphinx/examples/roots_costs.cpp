// The program of docs/sphinx/roots/costs.md: the evaluations per method in the page's tables. It asserts
// each count, typed here from the independent prototypes (docs/roots/design.md, Evidence, Phases 1 to 4 and 6, and
// docs/roots/PHASE6_API.md; Bisection on 1/x - 1 derived: half the width 2.5 * 2^-(k+1) against tol near 1, about
// 2^-51, gives k = 52, so 54, which .scratch/roots1/proto3.py also gives), each value within abs_error + c of its
// closed-form root, and abs_error <= 2 tol(value). The counts hold for fn written as here, in double, without FMA
// contraction. The page shows no code from this program.

#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <array>
#include <cmath>
#include <expected>
#include <limits>
#include <print>
#include <tuple>

namespace {

    namespace rt = nxx::roots;

    // Each case as fn returning (f, f', f''), in the prototypes' expressions and order; the methods without
    // derivatives use the first element.
    constexpr auto square = [](double x) { return std::tuple{x * x - 2.0, 2.0 * x, 2.0}; };
    constexpr auto wallis = [](double x) { return std::tuple{(x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0, 6.0 * x}; };
    constexpr auto x20 = [](double x) {
        double const x2 = x * x;
        double const x4 = x2 * x2;
        double const x8 = x4 * x4;
        double const x16 = x8 * x8;
        return std::tuple{x8 * x8 * x4 - 1.0, 20.0 * x16 * x2 * x, 380.0 * x16 * x2};
    };
    constexpr auto triple = [](double x) {
        return std::tuple{(x - 1.0) * (x - 1.0) * (x - 1.0), 3.0 * (x - 1.0) * (x - 1.0), 6.0 * (x - 1.0)};
    };
    constexpr auto reciprocal = [](double x) { return std::tuple{1.0 / x - 1.0, -1.0 / (x * x), 2.0 / (x * x * x)}; };

    int failures = 0;

    // fn's calls by solve<M> from start (a bracket or a guess); the result within abs_error + c of root and
    // abs_error <= 2 tol(value), else a failure.
    template<typename M, typename F, typename S>
    int calls_of(F const &fn, S const &start, double root, double c)
    {
        int calls = 0;
        auto const counted = [&](double x) {
            ++calls;
            return fn(x);
        };
        auto const r = rt::solve<M>(counted, start);
        double const eps = std::numeric_limits<double>::epsilon();
        if (!(r && std::abs(r->value - root) <= r->abs_error + c
              && r->abs_error <= 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * std::abs(r->value)))) {
            std::println("inaccurate: root {}", root);
            ++failures;
        }
        return calls;
    }

    // The calls of a value-only fn by Newton with nxx::deriv::derivative_of (Central<2>), with the same checks.
    template<typename F, typename S>
    int numerical_calls_of(F const &fn, S const &start, double root, double c)
    {
        int calls = 0;
        auto const value = [&](double x) {
            ++calls;
            return std::get<0>(fn(x));
        };
        auto const r = rt::solve<rt::Newton>(value, nxx::deriv::derivative_of(value), start);
        double const eps = std::numeric_limits<double>::epsilon();
        if (!(r && std::abs(r->value - root) <= r->abs_error + c
              && r->abs_error <= 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * std::abs(r->value)))) {
            std::println("inaccurate (numerical): root {}", root);
            ++failures;
        }
        return calls;
    }

    void expect(int got, int want, char const *what)
    {
        if (got != want) {
            std::println("{}: {} evaluations, the table says {}", what, got, want);
            ++failures;
        }
    }

    // One row of the bracket table: Bisection, Brent, Toms748, AndersonBjorck, Ridders, Illinois, RegulaFalsi,
    // Newton, Halley, Newton with derivative_of.
    template<typename F>
    void row(char const *name, F const &fn, rt::Bracket<double> b, double root, double c, std::array<int, 10> want)
    {
        std::array<int, 10> const got{calls_of<rt::Bisection>(fn, b, root, c), calls_of<rt::Brent>(fn, b, root, c),
                                      calls_of<rt::Toms748>(fn, b, root, c),
                                      calls_of<rt::AndersonBjorck>(fn, b, root, c),
                                      calls_of<rt::Ridders>(fn, b, root, c), calls_of<rt::Illinois>(fn, b, root, c),
                                      calls_of<rt::RegulaFalsi>(fn, b, root, c), calls_of<rt::Newton>(fn, b, root, c),
                                      calls_of<rt::Halley>(fn, b, root, c), numerical_calls_of(fn, b, root, c)};
        for (std::size_t i = 0; i < got.size(); ++i) expect(got[i], want[i], name);
    }

} // namespace

int main()
{
    // c per case, fn's rounding near the root over |f'| plus half an ulp of the root for the reference in a double:
    // x^2 - 2: 7.9e-17 + 1.1e-16; x^3 - 2x - 5: 2e-16 + 2.3e-16; x^20 - 1: 1.1e-16 + 1.1e-16; (x - 1)^3: 0 (x - 1
    // is exact near 1 and each product keeps its sign) + 0 (the root 1 is a double); 1/x - 1: 1.2e-16 + 1.1e-16.
    // A numerical derivative leaves the bound alone (the bound never depends on the derivative).
    row("x^2 - 2", square, {1.0, 2.0}, 1.4142135623730950488, 2e-16, {52, 9, 9, 9, 12, 11, 20, 8, 6, 24});
    row("x^3 - 2x - 5", wallis, {2.0, 3.0}, 2.0945514815423265915, 4.3e-16, {51, 8, 9, 9, 14, 10, 23, 7, 6, 21});
    row("x^20 - 1", x20, {0.0, 5.0}, 1.0, 2.2e-16, {55, 19, 21, 30, 18, 30, 45, 11, 11, 33});
    row("(x - 1)^3", triple, {0.0, 3.0}, 1.0, 0.0, {54, 156, 149, 152, 86, 148, 145, 91, 54, 549});
    row("1/x - 1", reciprocal, {0.5, 3.0}, 1.0, 2.3e-16, {54, 11, 6, 5, 12, 12, 25, 9, 6, 27});

    // The guess table: x^2 - 2 from 1, by Newton, Halley, Secant, Steffensen, and Newton with derivative_of.
    double const s2 = 1.4142135623730950488;
    expect(calls_of<rt::Newton>(square, 1.0, s2, 2e-16), 7, "Newton from 1");
    expect(calls_of<rt::Halley>(square, 1.0, s2, 2e-16), 5, "Halley from 1");
    expect(calls_of<rt::Secant>(square, 1.0, s2, 2e-16), 9, "Secant from 1");
    expect(calls_of<rt::Steffensen>(square, 1.0, s2, 2e-16), 10, "Steffensen from 1");
    expect(numerical_calls_of(square, 1.0, s2, 2e-16), 21, "Newton with derivative_of from 1");
    // The footnote: on 1000 (x^2 - 2) from 1, Steffensen finds -sqrt(2) (docs/roots/PHASE6_API.md: 27 evaluations).
    // fn's rounding is 1000 times x^2 - 2's and so is |f'|, so c is the same as for x^2 - 2.
    auto const scaled = [](double x) { return 1000.0 * (x * x - 2.0); };
    expect(calls_of<rt::Steffensen>(scaled, 1.0, -s2, 2e-16), 27, "Steffensen on 1000 (x^2 - 2) from 1");
    return failures == 0 ? 0 : 1;
}
