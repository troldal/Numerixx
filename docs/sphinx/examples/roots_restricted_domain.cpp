// The program of docs/sphinx/how-to/roots-restricted-domain.md: a fn defined on part of the line. The page shows the
// parts between the markers; main returns non-zero unless what the page says holds, against closed forms, mpmath, or
// the independent prototype's search (.scratch/roots5/proto5d.py, as in tests/test_roots_search.cpp).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <cmath>
#include <limits>
#include <print>
#include <utility>
// </includes>
// clang-format on

#include <numbers>

int main()
{
    namespace rt = nxx::roots;
    double const eps = std::numeric_limits<double>::epsilon();

    // clang-format off
    // <leave>
    // ln x - 1, defined for x > 0; std::log gives NaN below 0. Its root is e.
    auto const fn = [](double x) { return std::log(x) - 1.0; };
    auto const fn_and_slope = [](double x) { return std::pair{std::log(x) - 1.0, 1.0 / x}; };
    auto const newton = rt::solve<rt::Newton>(fn_and_slope, 20.0);
    if (!newton) std::println("Newton from 20: {}", newton.error());
    auto const brent = rt::search(fn, 20.0).and_then([&](rt::Bracket<double> const &b) {
        return rt::solve<rt::Brent>(fn, b);
    });
    if (brent) std::println("search and Brent from 20: {:.17g}", brent->value);
    // </leave>

    // <nan>
    // Van der Waals' pressure of CO2 at 280 K, less 4 MPa, for the molar volume v in m^3/mol. The formula is finite
    // below the co-volume b = 4.27e-5, where it means nothing, so return NaN there.
    auto const pressure = [](double v) {
        if (!(v > 4.27e-5)) return std::numeric_limits<double>::quiet_NaN();
        return 8.314 * 280.0 / (v - 4.27e-5) - 0.364 / (v * v) - 4e6;
    };
    auto const volume = rt::search(pressure, 1e-4).and_then([&](rt::Bracket<double> const &b) {
        return rt::solve<rt::Brent>(pressure, b);
    });
    if (volume) std::println("v = {:.6e} m^3/mol", volume->value);
    // </nan>

    // <stencil>
    // ln(x / 1e-12), whose root is 1e-12: from the guess 1e-13 the default step of the numerical derivative, about
    // 6e-6, reaches below 0, so the derivative fails there; a bracket still converges.
    auto const g = [](double x) { return std::log(x / 1e-12); };
    auto const from_guess = rt::solve<rt::Newton>(g, nxx::deriv::derivative_of(g), 1e-13);
    auto const in_bracket = rt::solve<rt::Newton>(g, nxx::deriv::derivative_of(g), rt::Bracket{1e-20, 1.0});
    if (!from_guess) std::println("Newton from 1e-13: {}", from_guess.error());
    if (in_bracket) std::println("Newton in [1e-20, 1]: {:.17g}", in_bracket->value);
    // </stencil>
    // clang-format on

    // Newton's first step from 20 is -(ln 20 - 1) / (1/20) = -39.9, to x = -19.9, where fn is NaN.
    if (!(newton == std::unexpected(rt::error::non_finite))) return 1;
    // e within abs_error + c: ln x near e errs by 2 ulps of 1 (4.4e-16), - 1 is exact, over f'(e) = 1/e = 0.37:
    // 1.2e-15, plus half an ulp of e for the reference in a double, 2.2e-16.
    if (!brent || !(std::abs(brent->value - std::numbers::e) <= brent->abs_error + 1.2e-15 + 2.3e-16)) return 2;

    // The search from 1e-4: the prototype's 29 evaluations and bracket.
    int calls = 0;
    auto const b = rt::search(
        [&](double v) {
            ++calls;
            return pressure(v);
        },
        1e-4);
    if (!(b == rt::Bracket{0x1.78f49491f2dc3p-12, 0x1.00808c8259e20p-11} && calls == 29)) return 3;
    // v = 4.36658136461612262e-4 (mpmath, with the constants as the doubles above). Each term of fn errs by a few
    // units of its magnitude (5.9e6, 1.9e6, 4e6), about 8e-9 in all, over |f'| = 6.3e9: c = 1.3e-18; with half an ulp
    // of v, 2e-18.
    double const v = 4.36658136461612262e-4;
    if (!volume || !(std::abs(volume->value - v) <= volume->abs_error + 2e-18)) return 4;

    // The stencil: not_converged from the guess. In the bracket, the root 1e-12 within abs_error + c: x / 1e-12 near
    // 1 rounds by half an ulp of 1 and ln adds 2 ulps of its result, 5.5e-16 in all, over f' = 1e12: 5.5e-28; with
    // half an ulp of 1e-12, 7e-28.
    if (!(from_guess == std::unexpected(rt::error::not_converged))) return 5;
    if (!in_bracket || !(std::abs(in_bracket->value - 1e-12) <= in_bracket->abs_error + 7e-28)) return 6;
    if (!(in_bracket->abs_error <= 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * in_bracket->value)))
        return 7;
    return 0;
}
