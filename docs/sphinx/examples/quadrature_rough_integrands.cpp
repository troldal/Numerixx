// The program of docs/sphinx/quadrature/rough-integrands.md: singularities at an end, at an end other than 0 and
// inside the interval, a kink, and a feature that no method sees. The page shows the parts between the markers and the
// output in the raw strings; main returns non-zero unless the program prints that output and what the page says holds,
// against closed forms. The quoted integrands use only + - * /, std::sqrt, std::abs and comparisons, which IEEE 754
// makes exact or correctly rounded, so the output is the same on every target that does not contract a * b + c into a
// fused multiply-add; the two that call std::log and std::pow are described on the page by the bounds checked here.

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <limits>
#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
    namespace qd = nxx::quadrature;

    // clang-format off
    // <endpoint>
    int evaluations = 0;
    auto const inv_sqrt = [&evaluations](double x) {
        ++evaluations;
        return 1.0 / std::sqrt(x);
    };
    auto const k = qd::integrate(inv_sqrt, 0.0, 1.0); // the default, GaussKronrod<15>
    if (k) std::println("default:  {:.16f}, abs_error {:.1e}, {} evaluations", k->value, k->abs_error, evaluations);
    evaluations = 0;
    auto const t = qd::integrate<qd::TanhSinh>(inv_sqrt, 0.0, 1.0);
    if (t) std::println("TanhSinh: {:.16f}, abs_error {:.1e}, {} evaluations", t->value, t->abs_error, evaluations);
    // </endpoint>

    // <other_end>
    // 1 / sqrt(1 - x) over [0, 1], as written, and in u = 1 - x
    auto const in_x = qd::integrate<qd::TanhSinh>([](double x) { return 1.0 / std::sqrt(1.0 - x); }, 0.0, 1.0);
    auto const in_u = qd::integrate<qd::TanhSinh>([](double u) { return 1.0 / std::sqrt(u); }, 0.0, 1.0);
    if (in_x && in_u) std::println("in x: error {:.1e};  in u: error {:.1e}", std::abs(in_x->value - 2.0),
                                   std::abs(in_u->value - 2.0));
    // </other_end>

    // <kink>
    double const c = 1.0 / 3.0; // where the kink is
    auto const kink = [c](double x) { return std::abs(x - c); };
    auto const whole = qd::integrate(kink, 0.0, 1.0);
    auto const left = qd::integrate(kink, 0.0, c);
    auto const right = qd::integrate(kink, c, 1.0);
    if (whole) std::println("whole: {:.17f}, abs_error {:.1e}", whole->value, whole->abs_error);
    if (left && right) std::println("split: {:.17f}, abs_error {:.1e}", left->value + right->value,
                                    left->abs_error + right->abs_error);
    // </kink>

    // <inside>
    // 1 / sqrt|x - c| over [0, 1]: at once, then split at c with each half written in its distance u from c
    double const exact = 2.0 * std::sqrt(c) + 2.0 * std::sqrt(1.0 - c);
    auto const at_once = qd::integrate([c](double x) { return 1.0 / std::sqrt(std::abs(x - c)); }, 0.0, 1.0);
    if (at_once) std::println("at once: error {:.1e}, abs_error {:.1e}", std::abs(at_once->value - exact),
                              at_once->abs_error);
    auto const half = [](double u) { return 1.0 / std::sqrt(u); };
    auto const below = qd::integrate<qd::TanhSinh>(half, 0.0, c);       // x = c - u
    auto const above = qd::integrate<qd::TanhSinh>(half, 0.0, 1.0 - c); // x = c + u
    if (below && above) std::println("split:   error {:.1e}", std::abs(below->value + above->value - exact));
    // </inside>

    // <narrow>
    // a box of width 2e-4 around 0.51: its integral is 2e-4
    auto const box = [](double x) { return std::abs(x - 0.51) < 1e-4 ? 1.0 : 0.0; };
    auto const unseen = qd::integrate(box, 0.0, 1.0);
    if (unseen) std::println("value {}, abs_error {:.1e}", unseen->value, unseen->abs_error);
    // </narrow>
    // clang-format on

    constexpr std::string_view endpoint_output = R"(
default:  1.9999999862243827, abs_error 2.1e-08, 1350 evaluations
TanhSinh: 2.0000000000000000, abs_error 2.2e-14, 131 evaluations
)";
    constexpr std::string_view other_end_output = R"(
in x: error 1.0e-08;  in u: error 0.0e+00
)";
    constexpr std::string_view kink_output = R"(
whole: 0.27777777943220033, abs_error 2.4e-09
split: 0.27777777777777779, abs_error 3.1e-15
)";
    constexpr std::string_view inside_output = R"(
at once: error 5.4e-08, abs_error 4.1e-08
split:   error 0.0e+00
)";
    constexpr std::string_view narrow_output = R"(
value 0, abs_error 2.5e-321
)";
    std::string const expected = std::string(docs::quoted(endpoint_output))
                                 + std::string(docs::quoted(other_end_output)) + std::string(docs::quoted(kink_output))
                                 + std::string(docs::quoted(inside_output)) + std::string(docs::quoted(narrow_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!k || !t || !in_x || !in_u || !whole || !left || !right || !at_once || !below || !above || !unseen) return 2;
    double const eps = std::numeric_limits<double>::epsilon();

    // The endpoint singularity: both estimates cover their errors; "ten times fewer evaluations" (131 against 1,350).
    if (!(std::abs(k->value - 2.0) <= k->abs_error && std::abs(t->value - 2.0) <= t->abs_error)) return 3;
    // log x over [0, 1], -1: "about 110 evaluations by TanhSinh, about 600 by the default", both within their
    // estimates (114 and 630 measured with GCC 16; std::log need not be correctly rounded, hence the ranges).
    int ts_log = 0;
    int gk_log = 0;
    auto const tl = qd::integrate<qd::TanhSinh>([&ts_log](double x) { ++ts_log; return std::log(x); }, 0.0, 1.0);
    auto const gl = qd::integrate([&gk_log](double x) { ++gk_log; return std::log(x); }, 0.0, 1.0);
    if (!(tl && gl && std::abs(tl->value + 1.0) <= tl->abs_error && std::abs(gl->value + 1.0) <= gl->abs_error))
        return 4;
    if (!(ts_log >= 100 && ts_log <= 130 && gk_log >= 500 && gk_log <= 800)) return 5;

    // The other end: TanhSinh's points cannot come nearer to 1 than the spacing of doubles below it, eps / 2. The
    // sliver [1 - eps / 2, 1] holds 2 sqrt(eps / 2) = 2.1e-8 of the integral of (1 - x)^(-1/2), an estimate of what
    // written in x can lose, not a bound: the method's end term adds part of it back, and 1.0e-8 is measured. "Up to
    // about 2e-8": between a tenth of the sliver and the sliver. In u the singular end is 0, where doubles are dense.
    double const sliver = 2.0 * std::sqrt(eps / 2.0);
    if (!(std::abs(in_x->value - 2.0) > 0.1 * sliver && std::abs(in_x->value - 2.0) <= sliver)) return 6;
    if (!(std::abs(in_u->value - 2.0) <= 2.0 * eps)) return 7;
    // "Writing fn(1.0 - u) would not help": 1 - u rounds back onto the numbers near 1, so the error stays.
    auto const rounded =
        qd::integrate<qd::TanhSinh>([](double u) { return 1.0 / std::sqrt(1.0 - (1.0 - u)); }, 0.0, 1.0);
    if (!(rounded && std::abs(rounded->value - 2.0) > 0.1 * sliver)) return 16;

    // The default's limit: x^(-3/4), 4, which it "accepts with an error above both its estimate and the request"
    // (1.5e-8 relative): a false success; TanhSinh is within its estimate.
    auto const strong = qd::integrate([](double x) { return std::pow(x, -0.75); }, 0.0, 1.0);
    auto const strong_ts = qd::integrate<qd::TanhSinh>([](double x) { return std::pow(x, -0.75); }, 0.0, 1.0);
    if (!(strong && std::abs(strong->value - 4.0) > strong->abs_error && std::abs(strong->value - 4.0) > 1.5e-8 * 4.0))
        return 8;
    if (!(strong_ts && std::abs(strong_ts->value - 4.0) <= strong_ts->abs_error)) return 9;

    // The kink: the reference c^2 / 2 + (1 - c)^2 / 2 for the double c, rounded within 2 ulps. "The default copes, at
    // 330 evaluations against 150 for each half"; "TanhSinh fails honestly": not_converged, with a best.
    double const kink_exact = c * c / 2.0 + (1.0 - c) * (1.0 - c) / 2.0;
    if (!(std::abs(whole->value - kink_exact) <= whole->abs_error)) return 10;
    if (!(std::abs(left->value + right->value - kink_exact) <= 4.0 * eps * kink_exact)) return 11;
    // "8 digits" at once: a relative error between 1e-9 and 1e-8 (6e-9 measured).
    double const whole_relative = std::abs(whole->value - kink_exact) / kink_exact;
    if (!(whole_relative > 1e-9 && whole_relative < 1e-8)) return 17;
    int whole_count = 0;
    int left_count = 0;
    int right_count = 0;
    (void)qd::integrate([&](double x) { ++whole_count; return kink(x); }, 0.0, 1.0);
    (void)qd::integrate([&](double x) { ++left_count; return kink(x); }, 0.0, c);
    (void)qd::integrate([&](double x) { ++right_count; return kink(x); }, c, 1.0);
    if (!(whole_count == 330 && left_count == 150 && right_count == 150)) return 12;
    // "A kink exactly on a boundary of the first 10 panels costs nothing extra": |x^2 - 1| on [0, 2], kink at 1,
    // 150 evaluations, 2 within 2 ulps.
    int boundary_count = 0;
    auto const on_boundary = qd::integrate([&](double x) { ++boundary_count; return std::abs(x * x - 1.0); }, 0.0, 2.0);
    if (!(on_boundary && boundary_count == 150 && std::abs(on_boundary->value - 2.0) <= 4.0 * eps)) return 19;
    // "At a tight request an unsplit kink's estimate can fall short of its error": |x^2 - 2| on [0, 2] at rel_tol
    // 1e-13, 8 sqrt(2) / 3 - 4 / 3 (rounded within 4 ulps of 2.44), ok with abs_error 6.1e-14 below the actual error
    // 7.1e-14; split at std::sqrt(2.0), whose distance to sqrt 2 (below 1.6e-16) times the jump in the slope (4 sqrt 2)
    // bounds what the split leaves, under 1e-15: within the reference's 4 ulps, "the error is 0" as printed by %.1e.
    auto const kink2 = [](double x) { return std::abs(x * x - 2.0); };
    double const kink2_exact = 8.0 * std::sqrt(2.0) / 3.0 - 4.0 / 3.0;
    qd::Criteria<double> const tight{ .rel_tol = 1e-13 };
    auto const unsplit = qd::integrate(kink2, 0.0, 2.0, tight);
    if (!(unsplit && unsplit->abs_error < std::abs(unsplit->value - kink2_exact)
          && unsplit->abs_error >= 6.05e-14 && unsplit->abs_error < 6.15e-14
          && std::abs(unsplit->value - kink2_exact) >= 7.05e-14 && std::abs(unsplit->value - kink2_exact) < 7.15e-14))
        return 20;
    double const root2 = std::sqrt(2.0);
    auto const below2 = qd::integrate(kink2, 0.0, root2, tight);
    auto const above2 = qd::integrate(kink2, root2, 2.0, tight);
    if (!(below2 && above2 && std::abs(below2->value + above2->value - kink2_exact) <= 4.0 * eps * kink2_exact))
        return 21;
    auto const kink_ts = qd::integrate<qd::TanhSinh>(kink, 0.0, 1.0);
    if (!(!kink_ts && kink_ts.error().reason == qd::error::not_converged && kink_ts.error().best)) return 13;

    // The interior singularity at once: "an error above its estimate and above the request", 1.5e-8 relative.
    if (!(std::abs(at_once->value - exact) > at_once->abs_error && std::abs(at_once->value - exact) > 1.5e-8 * exact))
        return 14;

    // The box: both methods report success with 0, the integral 2e-4 unseen; TanhSinh too.
    // "None of the first 150 points hits it": the default stops on its mesh.
    auto const unseen_ts = qd::integrate<qd::TanhSinh>(box, 0.0, 1.0);
    if (!(unseen->value == 0.0 && unseen_ts && unseen_ts->value == 0.0)) return 15;
    int box_count = 0;
    (void)qd::integrate([&](double x) { ++box_count; return box(x); }, 0.0, 1.0);
    if (!(box_count == 150)) return 18;
    return 0;
}
