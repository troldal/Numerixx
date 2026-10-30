// The program of docs/sphinx/basics.md, its section on results: testing a std::expected before reading it, printing
// its error, and composing steps with and_then, transform, or_else and value_or. The page shows the parts between the
// markers and the output in the raw string; main returns non-zero unless the program prints that output and what the
// page says holds, against closed forms.

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <limits>
#include <print>
#include <utility>
// </includes>
// clang-format on

#include "printed.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <expected>
#include <string>
#include <type_traits>

// A derivative's Estimate converts to its value, so `if (*d)` compiles and tests the value against 0.
static_assert(std::is_convertible_v<nxx::deriv::Estimate<double>, double>);

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
    namespace rt = nxx::roots;

    // clang-format off
    // <test>
    auto const fn = [](double x) { return x * x - 2.0; };
    for (auto const bracket : {rt::Bracket{1.0, 2.0}, rt::Bracket{2.0, 3.0}}) {
        auto const r = rt::solve<rt::Brent>(fn, bracket);
        if (r) std::println("root {:.6f}", r->value);
        else std::println("no root: {}", r.error());
    }
    // </test>

    // <compose>
    // From a guess: search for a bracket, solve in it, keep the value, and fall back on NaN if a step fails.
    double const x = rt::search(fn, 1.5)
                         .and_then([&](rt::Bracket<double> const &b) { return rt::solve<rt::Brent>(fn, b); })
                         .transform([](rt::Root<double> const &root) { return root.value; })
                         .value_or(std::numeric_limits<double>::quiet_NaN());
    std::println("x = {:.6f}", x);
    // </compose>

    // <or-else>
    // Newton from a guess; if it fails, Brent in a bracket.
    auto const fdf = [](double t) { return std::pair{t * t - 2.0, 2.0 * t}; };
    auto const s = rt::solve<rt::Newton>(fdf, 0.0).or_else([&](rt::error) {
        return rt::solve<rt::Brent>(fn, rt::Bracket{0.0, 2.0});
    });
    if (s) std::println("s = {:.6f}", s->value);
    // </or-else>
    // clang-format on

    constexpr std::string_view test_output = R"(
root 1.414214
no root: no_sign_change
)";
    constexpr std::string_view compose_output = R"(
x = 1.414214
s = 1.414214
)";
    if (!(docs::printed() == std::string(docs::quoted(test_output)) + std::string(docs::quoted(compose_output))))
        return 1;

    // x and s are sqrt(2) within abs_error + c (tests/roots_check.hpp's bound for x^2 - 2: c = 8e-17 for x * x
    // rounding over |f'| = 2.83, plus 1.2e-16 for the reference in a double). Newton from 0 fails at once (f'(0) = 0),
    // which is why the fallback ran.
    double const sqrt2 = 1.4142135623730950488;
    auto const r = rt::search(fn, 1.5).and_then([&](rt::Bracket<double> const &b) {
        return rt::solve<rt::Brent>(fn, b);
    });
    if (!r || !(x == r->value) || !(std::abs(x - sqrt2) <= r->abs_error + 8e-17 + 1.2e-16)) return 2;
    if (!(rt::solve<rt::Newton>(fdf, 0.0) == std::unexpected(rt::error::not_converged))) return 3;
    if (!s || !(std::abs(s->value - sqrt2) <= s->abs_error + 8e-17 + 1.2e-16)) return 4;
    return 0;
}
