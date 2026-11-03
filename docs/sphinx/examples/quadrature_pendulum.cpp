// The program of docs/sphinx/quadrature/tutorial.md: the period of a pendulum at large amplitude, through the complete
// elliptic integral K, computed in two forms and checked against the arithmetic-geometric mean. The page shows the
// parts between the markers and the output in the raw strings; main returns non-zero unless the program prints that
// output and what the page says holds. The smooth form calls std::cos, which need not be correctly rounded, so its
// results are described on the page by the bounds checked here, and only its ratios, to six digits, are quoted; the
// algebraic form, the reference and the divergent case use only + - * / and std::sqrt, which IEEE 754 rounds correctly,
// so their output is the same on every target that does not contract a * b + c into a fused multiply-add.

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

#include <array>
#include <string>

namespace {

    namespace qd = nxx::quadrature;

    // clang-format off
    // <reference>
    // K = pi / (2 AGM(1, k')): a and b move to their arithmetic and geometric means until they agree to a few ulps;
    // the last step's mean, (a + b) / 2, is folded into the quotient
    double reference(double kp)
    {
        double const eps = std::numeric_limits<double>::epsilon();
        double a = 1.0;
        double b = kp;
        while (a - b > 4.0 * eps * a) {
            double const mean = (a + b) / 2.0;
            b = std::sqrt(a * b);
            a = mean;
        }
        return std::numbers::pi / (a + b);
    }
    // </reference>

    // <smooth>
    // K in its smooth form: the integral of 1 / sqrt(k'^2 + k^2 cos^2 phi) over [0, pi/2]
    auto smooth(double kp)
    {
        double const k2 = 1.0 - kp * kp; // k^2
        auto const fn = [kp, k2](double phi) {
            double const c = std::cos(phi);
            return 1.0 / std::sqrt(kp * kp + k2 * c * c);
        };
        return qd::integrate(fn, 0.0, std::numbers::pi / 2.0);
    }
    // </smooth>

    // <algebraic>
    // K in its algebraic form, with u = 1 - sin phi: the integral of 1 / sqrt(w (k'^2 + k^2 w)) over [0, 1],
    // w = u (2 - u), singular at u = 0
    auto algebraic(double kp, int &evaluations)
    {
        double const k2 = 1.0 - kp * kp;
        auto const fn = [kp, k2, &evaluations](double u) {
            ++evaluations;
            double const w = u * (2.0 - u);
            return 1.0 / (std::sqrt(w) * std::sqrt(kp * kp + k2 * w));
        };
        return qd::integrate<qd::TanhSinh>(fn, 0.0, 1.0);
    }
    // </algebraic>
    // clang-format on

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <first>
    double const kp = std::sqrt(0.5); // k' = cos(theta0 / 2) for theta0 = 90 degrees
    auto const K = smooth(kp);
    if (K) std::println("K = {:.15f}, abs_error {:.1e}", K->value, K->abs_error);
    // </first>

    // <check>
    double const exact = reference(kp);
    if (K) std::println("error {:.1e}", std::abs(K->value - exact));
    // </check>

    // <periods>
    for (double const degrees : {60.0, 90.0, 120.0}) {
        // k' = cos(theta0 / 2): cos 30, cos 45 and cos 60 degrees, written exactly
        double const c = degrees == 60.0 ? std::sqrt(3.0) / 2.0 : degrees == 90.0 ? std::sqrt(0.5) : 0.5;
        auto const k = smooth(c);
        if (k) std::println("{:5.0f} degrees: T / T0 = {:.6f}", degrees, 2.0 * k->value / std::numbers::pi);
    }
    // </periods>

    // <near>
    for (double const c : {1e-2, 1e-4, 1e-6, 1e-8}) {
        int evaluations = 0;
        auto const a = algebraic(c, evaluations);
        if (a) std::println("k' = {:.0e}: K = {:.15f}, error {:.1e}, {} evaluations", c, a->value,
                            std::abs(a->value - reference(c)), evaluations);
    }
    // </near>

    // <divergent>
    int evaluations = 0;
    auto const top = algebraic(0.0, evaluations);
    if (!top) std::println("{}, best {}", top.error(), top.error().best ? "present" : "none");
    // </divergent>
    // clang-format on

    constexpr std::string_view periods_output = R"(
   60 degrees: T / T0 = 1.073182
   90 degrees: T / T0 = 1.180341
  120 degrees: T / T0 = 1.372881
)";
    constexpr std::string_view near_output = R"(
k' = 1e-02: K = 5.991589340506995, error 8.9e-16, 261 evaluations
k' = 1e-04: K = 10.596634757087664, error 1.8e-15, 265 evaluations
k' = 1e-06: K = 15.201804919087710, error 5.3e-15, 526 evaluations
k' = 1e-08: K = 19.806975105072254, error 0.0e+00, 531 evaluations
)";
    constexpr std::string_view divergent_output = R"(
not_converged, best none
)";
    // The first two lines go through std::cos: checked below by bounds, and skipped in the comparison.
    std::string const all = docs::printed();
    std::size_t const third = all.find('\n', all.find('\n') + 1) + 1;
    std::string const expected = std::string(docs::quoted(periods_output)) + std::string(docs::quoted(near_output))
                                 + std::string(docs::quoted(divergent_output));
    if (!(all.substr(third) == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!K || top) return 2;
    double const eps = std::numeric_limits<double>::epsilon();
    // The reference's own error (numerics reviews): each step's rounding moves the AGM's limit by at most 1.5 u
    // relative, u = eps / 2, since the AGM is homogeneous of degree 1 and increasing in both arguments; the first
    // step's product 1 k' is exact, so its geometric mean costs only the square root's u / 2. The final quotient costs
    // about 2.35 u (fl(pi) 0.35 u, the sum and the division u each). With n = 8 steps at k' = 1e-8, the most on this
    // page, that is at most 8 (1.5 u) - u + 2.35 u = 13.35 u < 7 eps.
    double const agm_bound = 7.0 * eps;
    // "K = 1.854074677301372"; "the error is within abs_error"; "a few units in the last place at most": within the
    // reference's bound twice over, 14 eps relative; "the method met it with its first 150 evaluations", and
    // abs_error at its floor, about 2e-14.
    if (!(std::abs(K->value - exact) <= K->abs_error + agm_bound * exact)) return 3;
    if (!(std::abs(K->value - 1.854074677301372) <= 1e-15)) return 4;
    if (!(K->abs_error < 3e-14 && std::abs(K->value - exact) <= 2.0 * agm_bound * exact)) return 5;
    int first = 0;
    auto const counted_first = qd::integrate([&first](double phi) {
        ++first;
        double const c = std::cos(phi);
        return 1.0 / std::sqrt(0.5 + 0.5 * c * c);
    }, 0.0, std::numbers::pi / 2.0);
    if (!(counted_first && first == 150)) return 13;
    // The smooth form near 180 degrees. pi / 2 rounds to a double delta = 6.123233995736766e-17 below it, so the
    // interval misses [fl(pi / 2), pi / 2], where the integrand is 1 / sqrt(k'^2 + k^2 sin^2 psi) ~ 1 / k' (psi < delta
    // << k'): a sliver of delta / k' to within a relative 1e-8. The method's points near fl(pi / 2) are rounded by up
    // to half an ulp, 1.1e-16, where the integrand changes by up to 1 / k' over the panel: an error of the same order,
    // either sign (numerics review). "At k' = 1e-6 the result is off by more than 1e-11, at 1e-8 by more than 1e-9"
    // (8.1e-11 and 3.8e-9 measured); "at 1e-6 most of that is the sliver"; "both within the request and abs_error".
    double const delta = 6.123233995736766e-17;
    auto const s6 = smooth(1e-6);
    auto const s8 = smooth(1e-8);
    if (!(s6 && s8 && std::abs(s6->value - reference(1e-6)) > 1e-11 && std::abs(s8->value - reference(1e-8)) > 1e-9))
        return 6;
    double const after = std::abs(s6->value + delta / 1e-6 - reference(1e-6));
    if (!(after < 0.5 * std::abs(s6->value - reference(1e-6)) && after <= s6->abs_error)) return 7;
    if (!(std::abs(s6->value - reference(1e-6)) <= s6->abs_error
          && std::abs(s8->value - reference(1e-8)) <= s8->abs_error))
        return 14;
    // "Asked for 1e-14 at k' = 1e-4, the method reports an abs_error below the error against K, and yet above the
    // error against the integral over [0, fl(pi / 2)]": the estimate is about the interval it was given.
    double const k4 = 1.0 - 1e-8;
    auto const tight = qd::integrate([k4](double phi) {
        double const c = std::cos(phi);
        return 1.0 / std::sqrt(1e-8 + k4 * c * c);
    }, 0.0, std::numbers::pi / 2.0, {.rel_tol = 1e-14});
    double const K4 = reference(1e-4);
    if (!(tight && tight->abs_error < std::abs(tight->value - K4)
          && std::abs(tight->value - (K4 - delta / 1e-4)) <= tight->abs_error))
        return 11;
    // "Measure the angle from the top instead": psi = pi / 2 - phi, cos phi = sin psi, the peak at the exact end 0;
    // within 1e-13 of K at k' = 1e-6 and 1e-8 (0 and 2.8e-14 measured).
    for (double const c : {1e-6, 1e-8}) {
        double const k2c = 1.0 - c * c;
        auto const top_first = qd::integrate([c, k2c](double psi) {
            double const sn = std::sin(psi);
            return 1.0 / std::sqrt(c * c + k2c * sn * sn);
        }, 0.0, std::numbers::pi / 2.0);
        if (!(top_first && std::abs(top_first->value - reference(c)) <= 1e-13)) return 12;
    }
    // "The default needs four to eleven times as many evaluations on the algebraic form": measured 10.3 at 90 degrees
    // (1,350 against 131), 6.2, 7.4, 4.4 and 5.1 at the k' printed.
    for (double const c : {std::sqrt(0.5), 1e-2, 1e-4, 1e-6, 1e-8}) {
        int ts = 0;
        int gk = 0;
        auto const counted_ts = algebraic(c, ts);
        double const k2 = 1.0 - c * c;
        auto const counted_gk = qd::integrate([c, k2, &gk](double u) {
            ++gk;
            double const w = u * (2.0 - u);
            return 1.0 / (std::sqrt(w) * std::sqrt(c * c + k2 * w));
        }, 0.0, 1.0);
        if (!(counted_ts && counted_gk && gk >= 4 * ts && gk <= 11 * ts)) return 8;
    }
    // "Written in t, with the singularity at t = 1, TanhSinh is 8e-9 off at 60 degrees": the spacing of doubles near 1.
    double const k60 = 0.25;
    auto const naive = qd::integrate<qd::TanhSinh>([k60](double t) {
        return 1.0 / std::sqrt((1.0 - t * t) * (1.0 - k60 * t * t));
    }, 0.0, 1.0);
    double const ref60 = reference(std::sqrt(3.0) / 2.0);
    if (!(naive && std::abs(naive->value - ref60) > 5e-9 && std::abs(naive->value - ref60) < 1.5e-8)) return 9;
    // The divergent case: 5,159 evaluations, TanhSinh's ten levels.
    if (!(top.error().reason == qd::error::not_converged && !top.error().best && evaluations == 5159)) return 10;
    return 0;
}
