// The program of docs/sphinx/quadrature/method-of-your-own.md: a composite Simpson rule with a run-time number of
// panels, written as a fixed rule, checked against the concept where it is defined, and called through integrate() on
// 1 / (1 + x) over [0, 1], whose integral is ln 2. The page shows the parts between the markers and the output in the
// raw string; main returns non-zero unless the program prints that output and what the page says holds. The rule and
// the integrand use only + - * /, which IEEE 754 rounds correctly, so the output is the same on every target that does
// not contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <cassert>
#include <cmath>
#include <expected>
#include <functional>
#include <limits>
#include <numbers>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

// clang-format off
// <rule>
namespace qd = nxx::quadrature;

// The composite Simpson rule on `panels` equal panels of [a, b]: each panel integrated from its ends and its
// midpoint, with the weights h/6, 4h/6 and h/6, h the panel's width
struct Simpson
{
    int panels = 10;

    template<nxx::real T, typename F>
    constexpr std::expected<T, qd::error> apply(F &fn, T const &a, T const &b) const
    {
        assert(panels >= 1);      // a precondition: fewer panels is a bug in the caller, not an error to report
        T const width = T(b - a); // a < b, both finite, but b - a may overflow
        if (!nxx::math::isfinite(width)) return std::unexpected(qd::error::non_finite);
        T const h = T(width / T(panels));
        // The point a + t h, kept inside [a, b]: rounding could carry the last ones past b
        auto const at = [&](T const &t) {
            T const x = T(a + T(t * h));
            return b < x ? b : x;
        };
        auto const f = [&](T const &x) { return T(std::invoke(fn, x)); };
        T ends = T(f(a) + f(b)); // the outer ends, weight 1
        T inner = T(0);          // the ends shared by two panels, weight 2
        T middles = T(0);        // the panels' midpoints, weight 4
        for (int i = 0; i < panels; ++i) {
            if (i > 0) inner = T(inner + f(at(T(i))));
            middles = T(middles + f(at(T(T(i) + T(1) / T(2)))));
        }
        return T(T(T(ends + T(T(2) * inner)) + T(T(4) * middles)) * T(h / T(6)));
    }
};
static_assert(qd::fixed_rule<Simpson, double>); // checked where it is defined
// </rule>
// clang-format on

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <call>
    auto const fn = [](double x) { return 1.0 / (1.0 + x); };
    for (int const n : {10, 100}) {
        auto const s = qd::integrate(Simpson{.panels = n}, fn, 0.0, 1.0);
        if (s) std::println("{:3} panels: {:.17g}, error {:.1e}", n, *s, std::abs(*s - std::numbers::ln2));
    }
    // </call>
    // clang-format on

    constexpr std::string_view call_output = R"(
 10 panels: 0.69314737466511611, error 1.9e-07
100 panels: 0.69314718057947522, error 2.0e-11
)";
    if (!(docs::printed() == docs::quoted(call_output))) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    auto const s10 = qd::integrate(Simpson{.panels = 10}, fn, 0.0, 1.0);
    auto const s100 = qd::integrate(Simpson{.panels = 100}, fn, 0.0, 1.0);
    if (!s10 || !s100) return 2;
    // "The error falls ten thousand times for ten times as many panels". Simpson's rule is (4 T(h/2) - T(h)) / 3 for
    // the trapezoidal sums T, whose Euler-Maclaurin expansion T - I = h^2/12 D1 - h^4/720 D3 + h^6/30240 D5
    // - h^8/1209600 D7 + ... (Dk = f^(k)(b) - f^(k)(a)) gives S - I = h^4/2880 D3 - h^6/96768 D5 + 0.328 h^8/1209600 D7
    // - ... for panels of width h. With f = 1 / (1 + x), D3 = 5.625, D5 = 118.125 and D7 = 5020.3, so
    // S - I = 1.953125e-3 h^4 - 1.2207e-3 h^6 + 1.36e-3 h^8 - ...; the alternating tail lies below 2e-3 h^8 (1.3e-11
    // measured at h = 0.1, against 1.36e-11). The rounding: three plain sums of at most n terms of size at most 1,
    // weighted 1, 2 and 4 and scaled by h/6 = 1 / (6n), err by at most about n eps / 2 (the ln 2 constant adds half an
    // ulp); n eps bounds it. So each error lies within 2e-3 h^8 + n eps of the two-term expansion: 2.0e-11 at n = 10,
    // 2.2e-14 at n = 100, against errors of 1.9e-7 and 2.0e-11.
    double const eps = std::numeric_limits<double>::epsilon();
    auto const expansion = [](double h) { return 1.953125e-3 * h * h * h * h - 1.2207e-3 * h * h * h * h * h * h; };
    auto const within = [&](double error, double h, int n) {
        return std::abs(error - expansion(h)) <= 2e-3 * h * h * h * h * h * h * h * h + n * eps;
    };
    if (!(within(*s10 - std::numbers::ln2, 0.1, 10) && within(*s100 - std::numbers::ln2, 0.01, 100)))
        return 3;
    // integrate() puts the ends in order and validates them before the rule sees them: reversed ends give the negated
    // value, equal ends 0 without calling fn, an infinite end non_finite.
    if (!(qd::integrate(Simpson{.panels = 10}, fn, 1.0, 0.0) == -*s10)) return 4;
    int calls = 0;
    auto const counted = [&calls](double x) {
        ++calls;
        return x;
    };
    if (!(qd::integrate(Simpson{}, counted, 0.5, 0.5) == 0.0 && calls == 0)) return 5;
    double const inf = std::numeric_limits<double>::infinity();
    if (!(qd::integrate(Simpson{}, fn, 0.0, inf) == std::unexpected(qd::error::non_finite))) return 6;
    // "A value that is not finite becomes non_finite": 1/x has 1/0 = inf at the end a = 0, which Simpson evaluates.
    if (!(qd::integrate(Simpson{}, [](double x) { return 1.0 / x; }, 0.0, 1.0)
          == std::unexpected(qd::error::non_finite)))
        return 7;
    // "It runs at compile time too", since apply is constexpr and the integrand is: x^2 on [0, 3] is 9, exact for
    // Simpson's rule (degree 3), and every operation is exact in double for one panel: (0 + 9 + 4 * 2.25) * 0.5 = 9.
    static_assert(qd::integrate(Simpson{.panels = 1}, [](double x) { return x * x; }, 0.0, 3.0) == 9.0);
    // The note: integrate<Simpson>(fn, a, b) uses Simpson{}, whose default is 10 panels.
    if (!(qd::integrate<Simpson>(fn, 0.0, 1.0) == s10 && qd::integrate(Simpson{}, fn, 0.0, 1.0) == s10)) return 8;
    // The other concepts: Simpson is a fixed rule only.
    static_assert(!qd::iterative_method<Simpson, double> && !qd::estimating_method<Simpson, double>);
    return 0;
}
