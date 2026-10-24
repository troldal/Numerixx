// E18 (Chapter 7, section 7.5): is a stencil's result good enough when it comes without an error estimate? Four
// checks on the thread's cp at T0 = 1000 K: the same rule at h, 2h and h/2; a rule of higher order; Ridders as an
// independent estimate; and, where f' is near 0, the absolute rather than the relative error (sin x at pi/2, since
// the thread's derivatives have no zero in its range). The closed form is the answer key.

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <print>

using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Ridders;

namespace {

    double H(double T) { return book::shomate_H(book::nitrogen, T); }

} // namespace

// clang-format off
// <check-digits>
// The number of decimal digits on which a and b agree.
double digits(double a, double b)
{
    if (a == b) return 16.0;
    return -std::log10(std::abs(a - b) / std::abs(b));
}
// </check-digits>
// clang-format on

int main()
{
    bool ok = true;
    double const T0 = book::thread_T0;
    double const exact = book::shomate_cp(book::nitrogen, T0);

    // clang-format off
// <check-by-steps>
// Central<2> at its default step h, and at 2h and h/2
// the default step: root_epsilon<P + N, T>() * max(|x|, 1)
double const s_T = std::max(std::abs(T0), 1.0);
double const h = nxx::math::root_epsilon<3, double>() * s_T;
auto const d = derivative(H, T0, h);
auto const d2 = derivative(H, T0, 2 * h);
auto const d0 = derivative(H, T0, h / 2);
if (!d || !d2 || !d0) return 1;
double const agree = std::min(
    digits(d2->value, d->value),
    digits(d0->value, d->value));
std::println("h, 2h, h/2 agree to {:.1f} digits", agree);
// A physical check: a heat capacity is positive.
if (d->value <= 0) std::println("cp <= 0: wrong");
// </check-by-steps>
    // clang-format on
    std::println("  d = {:.12f}, 2h {:.12f}, h/2 {:.12f}, exact {:.12f}", d->value, d2->value, d0->value, exact);
    std::println("  agreement with the closed form: {:.1f} digits", digits(*d, exact));
    // The text claims at least 8 agreeing digits, among the three values and with the closed form: the model gives
    // about 4e-11 relative at h (eps^(2/3)), up to about 1e-10 at h/2, where rounding doubles (DS1).
    ok = ok && agree >= 8 && digits(*d, exact) >= 8 && digits(*d2, exact) >= 8 && digits(*d0, exact) >= 8;
    ok = ok && d->value > 0;

    // clang-format off
// <check-by-order>
// A rule of higher order, with its own default step
auto const c4 = derivative<Central<4>>(H, T0);
if (c4) {
    std::println("Central<4> agrees to {:.1f} digits",
        digits(d->value, c4->value));
}
// </check-by-order>
    // clang-format on
    // Central<4> is about 2e-13 off (DS1b), so the agreement measures the error of Central<2>: at least 8 digits.
    ok = ok && c4 && digits(*d, *c4) >= 8 && digits(*c4, exact) >= 8;

    // clang-format off
// <check-by-ridders>
// Ridders: an independent value with an error estimate.
// If that estimate is small, |d - r| is the error of d.
auto const r = derivative<Ridders>(H, T0);
if (r) {
    double const diff = std::abs(d->value - r->value);
    std::println("|d - r| = {:.1e}, Ridders +- {:.1e}",
        diff, *r->abs_error);
}
// </check-by-ridders>
    // clang-format on
    if (!r) return 1;
    double const diff = std::abs(d->value - r->value);
    double const actual = std::abs(d->value - exact);
    std::println("  the stencil's actual error {:.1e}; Ridders' actual error {:.1e}", actual,
        std::abs(r->value - exact));
    // Ridders estimates its own error, not the stencil's: its estimate (about 3e-12, DS1b) is far below the stencil's
    // error (about 8e-10). But when Ridders is covered (|r - exact| <= 6 * its estimate, the tests' criterion), the
    // triangle inequality makes |d - r| the stencil's error to within 6 times that estimate.
    ok = ok && std::abs(r->value - exact) <= 6 * *r->abs_error;
    ok = ok && std::abs(diff - actual) <= 6 * *r->abs_error;
    if (*r->abs_error < actual)
        std::println("  (Ridders' estimate is below the stencil's error: it estimates Ridders' own value)");

    // clang-format off
// <check-absolute>
// Near pi/2, f' = cos x is near 0: few digits of f' agree,
// but the spread is tiny against |f| / max(|x|, 1).
auto const sine = [](double x) { return std::sin(x); };
double const x = 1.5708; // just above pi/2
double const hs = nxx::math::root_epsilon<3, double>() * x;
auto const s = derivative(sine, x, hs);
auto const s2 = derivative(sine, x, 2 * hs);
auto const s0 = derivative(sine, x, hs / 2);
if (!s || !s2 || !s0) return 1;
double const spread = std::max(std::abs(*s2 - *s),
                               std::abs(*s0 - *s));
double const of_f1 = -std::log10(spread / std::abs(*s));
double const scale = std::abs(sine(x)) / x;
double const of_scale = -std::log10(spread / scale);
std::println("digits: {:.1f} of f', {:.1f} of the scale",
    of_f1, of_scale);
// </check-absolute>
    // clang-format on
    double const cos_x = std::cos(x);
    double const absolute = std::abs(s->value - cos_x);
    std::println("  f' = {:.6e}, cos x = {:.6e}; absolute error {:.1e}, relative {:.1e}", s->value, cos_x, absolute,
        absolute / std::abs(cos_x));
    // The absolute error is of the model's size, eps^(2/3) |f| / max(|x|, 1), about 4e-11: at most 1e-9 and at least
    // 8 digits of the scale. Relative to f' = -3.7e-6 the same error leaves far fewer: the text's point is
    // that the scale gives more agreeing digits than f' itself.
    ok = ok && absolute <= 1e-9 && of_scale >= 8 && of_f1 < of_scale;

    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}
