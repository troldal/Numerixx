// The program of docs/sphinx/quadrature/infinite-intervals.md: TanhSinh on the line and on half lines, the other
// methods' answer to an infinite end, a split between a finite and an infinite piece, and the limits of TanhSinh's
// maps. The page shows the parts between the markers and the output in the raw strings; main returns non-zero unless
// the program prints that output and what the page says holds, against closed forms. The quoted integrands use only
// + - * / and comparisons, which IEEE 754 rounds correctly, so the output is the same on every target that does not
// contract a * b + c into a fused multiply-add; those that call std::exp, std::log, std::pow and std::sin are described
// on the page by the bounds checked here.

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
    namespace qd = nxx::quadrature;
    double const pi = std::numbers::pi;

    // clang-format off
    // <line>
    double const inf = std::numeric_limits<double>::infinity();
    int evaluations = 0;
    auto const lorentz = [&evaluations](double x) {
        ++evaluations;
        return 1.0 / (1.0 + x * x);
    };
    auto const line = qd::integrate<qd::TanhSinh>(lorentz, -inf, inf); // pi
    if (line) std::println("line:      {:.16f}, abs_error {:.1e}, {} evaluations", line->value, line->abs_error,
                           evaluations);
    evaluations = 0;
    auto const half = qd::integrate<qd::TanhSinh>(lorentz, 0.0, inf); // pi / 2
    if (half) std::println("half line: {:.16f}, abs_error {:.1e}, {} evaluations", half->value, half->abs_error,
                           evaluations);
    // </line>

    // <others>
    auto const k = qd::integrate(lorentz, 0.0, inf); // the default, GaussKronrod<15>
    if (!k) std::println("no integral: {}", k.error());
    // </others>

    // <split>
    // 1 / (1 + x^2) on [0, 1), twice that beyond: 3 pi / 4
    auto const stepped = [](double x) { return (x < 1.0 ? 1.0 : 2.0) / (1.0 + x * x); };
    auto const at_once = qd::integrate<qd::TanhSinh>(stepped, 0.0, inf);
    if (!at_once) std::println("at once: {}", at_once.error());
    auto const near = qd::integrate(stepped, 0.0, 1.0);
    auto const far = qd::integrate<qd::TanhSinh>(stepped, 1.0, inf);
    if (near && far) std::println("split:   {:.16f}, abs_error {:.1e}", near->value + far->value,
                                  near->abs_error + far->abs_error);
    // </split>

    // <far_peak>
    // a peak at x0 = 40 on the line: at once, split at x0, and with u = x - x0
    double const x0 = 40.0;
    auto const peak = [x0](double x) { return 1.0 / (1.0 + (x - x0) * (x - x0)); };
    auto const whole = qd::integrate<qd::TanhSinh>(peak, -inf, inf);
    if (!whole) std::println("at once:       {}", whole.error());
    auto const below = qd::integrate<qd::TanhSinh>(peak, -inf, x0);
    auto const above = qd::integrate<qd::TanhSinh>(peak, x0, inf);
    if (below && above) std::println("split at x0:   error {:.1e}", std::abs(below->value + above->value - pi));
    auto const shifted = qd::integrate<qd::TanhSinh>([](double u) { return 1.0 / (1.0 + u * u); }, -inf, inf);
    if (shifted) std::println("with u = x-x0: error {:.1e}", std::abs(shifted->value - pi));
    // </far_peak>

    // <far_end>
    // 1 / (1 + (x - a)^2) on [a, inf), pi / 2, for a far from 0
    for (double const a : {1e10, 1e15}) {
        auto const r = qd::integrate<qd::TanhSinh>([a](double x) { return 1.0 / (1.0 + (x - a) * (x - a)); }, a, inf);
        if (r) std::println("a = {:.0e}: error {:.1e}, abs_error {:.1e}", a, std::abs(r->value - pi / 2.0),
                            r->abs_error);
    }
    // </far_end>

    // <overflow>
    // the integral of x^51 e^-x over [0, inf) is 51!
    auto const product = qd::integrate<qd::TanhSinh>([](double x) { return std::pow(x, 51) * std::exp(-x); },
                                                     0.0, inf);
    auto const logs = qd::integrate<qd::TanhSinh>([](double x) { return std::exp(51.0 * std::log(x) - x); },
                                                  0.0, inf);
    // </overflow>
    // clang-format on

    constexpr std::string_view line_output = R"(
line:      3.1415926535897931, abs_error 3.5e-14, 133 evaluations
half line: 1.5707963267948966, abs_error 1.7e-14, 133 evaluations
)";
    constexpr std::string_view others_output = R"(
no integral: non_finite
)";
    constexpr std::string_view split_output = R"(
at once: not_converged
split:   2.3561944901923448, abs_error 2.6e-14
)";
    constexpr std::string_view far_peak_output = R"(
at once:       not_converged
split at x0:   error 8.9e-16
with u = x-x0: error 0.0e+00
)";
    constexpr std::string_view far_end_output = R"(
a = 1e+10: error 3.5e-08, abs_error 1.1e-06
a = 1e+15: error 1.2e-03, abs_error 7.9e-02
)";
    std::string const expected = std::string(docs::quoted(line_output)) + std::string(docs::quoted(others_output))
                                 + std::string(docs::quoted(split_output)) + std::string(docs::quoted(far_peak_output))
                                 + std::string(docs::quoted(far_end_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!line || !half || k || at_once || !near || !far || whole || !below || !above || !shifted) return 2;
    double const eps = std::numeric_limits<double>::epsilon();

    // The line and the half line: within their estimates.
    if (!(std::abs(line->value - pi) <= line->abs_error && std::abs(half->value - pi / 2.0) <= half->abs_error))
        return 3;
    // "Every other method, Romberg and GaussLegendre<N> included, returns non_finite: it compiles, and fails at run
    // time".
    auto const romberg = qd::integrate<qd::Romberg>(lorentz, 0.0, inf);
    auto const gauss = qd::integrate<qd::GaussLegendre<10>>(lorentz, 0.0, inf);
    if (!(k.error().reason == qd::error::non_finite && !romberg && romberg.error().reason == qd::error::non_finite
          && !gauss && gauss.error() == qd::error::non_finite))
        return 4;
    // The split: 3 pi / 4 to within 4 ulps; at once, not_converged.
    if (!(std::abs(near->value + far->value - 3.0 * pi / 4.0) <= 4.0 * eps * pi
          && at_once.error().reason == qd::error::not_converged))
        return 5;
    // The far peak: "at x0 = 10 it still succeeds, at about 2,000 evaluations"; at 40 not_converged; split and
    // shifted within 4 ulps of pi.
    int ten = 0;
    auto const at_ten = qd::integrate<qd::TanhSinh>([&ten](double x) {
        ++ten;
        return 1.0 / (1.0 + (x - 10.0) * (x - 10.0));
    }, -inf, inf);
    if (!(at_ten && std::abs(at_ten->value - pi) <= at_ten->abs_error && ten > 1500 && ten < 2500)) return 6;
    if (!(whole.error().reason == qd::error::not_converged
          && std::abs(below->value + above->value - pi) <= 4.0 * eps * pi
          && std::abs(shifted->value - pi) <= 4.0 * eps * pi))
        return 7;
    // A half line from 1e17: "no point is inside" (a + 1 rounds to a), not_converged without a best; with u = x - a,
    // pi / 2 to 2 ulps.
    auto const beyond = qd::integrate<qd::TanhSinh>([](double x) { return 1.0 / (1.0 + (x - 1e17) * (x - 1e17)); },
                                                    1e17, inf);
    if (!(!beyond && beyond.error().reason == qd::error::not_converged && !beyond.error().best)) return 8;
    auto const moved = qd::integrate<qd::TanhSinh>([](double u) { return 1.0 / (1.0 + u * u); }, 0.0, inf);
    if (!(moved && std::abs(moved->value - pi / 2.0) <= 2.0 * eps * (pi / 2.0))) return 15;
    // Overflowing factors: e^-x is exactly 0 from x = 745 and x^51 inf from about 1.1e6, so beyond that the product
    // is inf * 0 = NaN. "The first points jump from about 230, where the integrand is still 2e20, to about 3.4e6": the
    // first value not finite is at x between 3e6 and 4e6, and the largest x sampled before it lies between 200 and
    // 300 with a value above 1e20. "With x^50 the same NaN is accepted": x^50 e^-x succeeds, within its estimate of
    // 50! (49 roundings in the product, 25 eps relative at most).
    double first_bad = 0.0;
    double last_x = 0.0;
    double last_value = 0.0;
    auto const watched = qd::integrate<qd::TanhSinh>([&](double x) {
        double const v = std::pow(x, 51) * std::exp(-x);
        if (first_bad == 0.0 && !std::isfinite(v)) first_bad = x;
        if (first_bad == 0.0 && x > last_x) {
            last_x = x;
            last_value = v;
        }
        return v;
    }, 0.0, inf);
    if (!(!watched && first_bad > 3e6 && first_bad < 4e6 && last_x > 200.0 && last_x < 300.0 && last_value > 1e20))
        return 16;
    double factorial50 = 1.0;
    for (int i = 2; i <= 50; ++i) factorial50 *= i;
    auto const fifty = qd::integrate<qd::TanhSinh>([](double x) { return std::pow(x, 50) * std::exp(-x); }, 0.0, inf);
    if (!(fifty && std::abs(fifty->value - factorial50) <= fifty->abs_error + 25.0 * eps * factorial50)) return 17;
    // Underflow is harmless after the terms are negligible: e^(-x^2/2) is 0 beyond 38.6 (x^2 / 2 > 745.1), and the
    // tail over [3, inf), sqrt(pi / 2) erfc(3 / sqrt 2) = 3.38369257395272758...e-3 (mpmath), is within its
    // estimate. An integrand tiny everywhere, 1e-300 x^-1.5 over [1, inf): not_converged.
    if (!(std::exp(-38.7 * 38.7 / 2.0) == 0.0 && std::exp(-38.5 * 38.5 / 2.0) > 0.0)) return 18;
    auto const tail = qd::integrate<qd::TanhSinh>([](double x) { return std::exp(-x * x / 2.0); }, 3.0, inf);
    if (!(tail && std::abs(tail->value - 3.3836925739527276e-3) <= tail->abs_error + 2.0 * eps * 3.4e-3)) return 19;
    auto const tiny = qd::integrate<qd::TanhSinh>([](double x) { return 1e-300 * std::pow(x, -1.5); }, 1.0, inf);
    if (!(!tiny && tiny.error().reason == qd::error::not_converged)) return 20;
    // Overflowing factors, as written on the page: the product is not finite where the tail rule
    // cannot accept it: non_finite. Written through logarithms, within its estimate of 51! (the product of the
    // integers, each multiplication exact or within half an ulp: 50 roundings, 50 eps / 2 relative at most).
    double factorial = 1.0;
    for (int i = 2; i <= 51; ++i) factorial *= i;
    if (!(!product && product.error().reason == qd::error::non_finite)) return 9;
    if (!(logs && std::abs(logs->value - factorial) <= logs->abs_error + 25.0 * eps * factorial)) return 10;
    // Slow decay: x^-1.01 over [1, inf), 100: "not_converged in double, with a best near 99.9 whose estimate covers
    // its error".
    auto const slow = qd::integrate<qd::TanhSinh>([](double x) { return std::pow(x, -1.01); }, 1.0, inf);
    if (!(!slow && slow.error().reason == qd::error::not_converged && slow.error().best
          && std::abs(slow.error().best->value - 100.0) < 0.2
          && std::abs(slow.error().best->value - 100.0) <= slow.error().best->abs_error))
        return 11;
    // A large number for infinity: [0, 1e300] not_converged, [0, inf) exact; [0, 1e10] is a different integral,
    // atan(1e10) = pi / 2 - 1e-10 to within 1e-30, which it gets within its estimate.
    auto const plain = [](double x) { return 1.0 / (1.0 + x * x); };
    auto const huge = qd::integrate<qd::TanhSinh>(plain, 0.0, 1e300);
    auto const wide = qd::integrate<qd::TanhSinh>(plain, 0.0, 1e10);
    if (!(!huge && huge.error().reason == qd::error::not_converged)) return 12;
    if (!(wide && std::abs(wide->value - (pi / 2.0 - 1e-10)) <= wide->abs_error)) return 13;
    // Oscillation: sin x / x over [0, inf), pi / 2: not_converged.
    auto const wave =
        qd::integrate<qd::TanhSinh>([](double x) { return x == 0.0 ? 1.0 : std::sin(x) / x; }, 0.0, inf);
    if (!(!wave && wave.error().reason == qd::error::not_converged)) return 14;
    return 0;
}
