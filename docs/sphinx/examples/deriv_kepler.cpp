// The program of docs/sphinx/deriv/tutorial.md: the derivatives of Kepler's eccentric anomaly E(M), a function computed
// by a root solve. The page compares them with the closed forms dE/dM = 1 / (1 - e cos E) and d2E/dM2 =
// -e sin E (dE/dM)^3 at the solved E; main checks them against constants computed independently of Numerixx (mpmath),
// and returns non-zero unless what the page says holds. The values go
// through std::sin, whose last digits differ between platforms, so the page quotes bounds, not output.

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <cmath>
#include <limits>
#include <print>

namespace dv = nxx::deriv;
namespace rt = nxx::roots;
// </includes>
// clang-format on

#include <format>
int main()
{
    // clang-format off
    // <function>
    double const e = 0.5; // the orbit's eccentricity
    // The eccentric anomaly E for a mean anomaly M: the root of E - e sin E - M in [M - e, M + e], or NaN if the
    // solve fails, which the derivative then reports as an error.
    auto const eccentric = [e](double M) {
        auto const kepler = [e, M](double E) { return E - e * std::sin(E) - M; };
        return rt::solve<rt::Brent>(kepler, rt::Bracket{M - e, M + e})
            .transform([](rt::Root<double> const &root) { return root.value; })
            .value_or(std::numeric_limits<double>::quiet_NaN());
    };
    // </function>

    // <central>
    double const M = 1.0;
    auto const slope = dv::derivative(eccentric, M);
    if (slope) std::println("dE/dM = {:.15f}", slope->value);
    // </central>

    // <exact>
    double const E = eccentric(M);
    double const exact = 1.0 / (1.0 - e * std::cos(E));
    if (slope) std::println("relative error {:.1e}", std::abs(slope->value - exact) / exact);
    // </exact>

    // <ridders>
    auto const r = dv::derivative<dv::Ridders>(eccentric, M);
    if (r) std::println("dE/dM = {:.15f}, abs_error {:.1e}, actual {:.1e}", r->value, *r->abs_error,
                        std::abs(r->value - exact));
    // </ridders>

    // <second>
    auto const r2 = dv::derivative<dv::Ridders, 2>(eccentric, M);
    double const exact2 = -e * std::sin(E) * exact * exact * exact;
    if (r2) std::println("d2E/dM2 = {:.12f}, abs_error {:.1e}, actual {:.1e}", r2->value, *r2->abs_error,
                         std::abs(r2->value - exact2));
    // </second>
    // clang-format on

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!slope || !r || !r->abs_error || !r2 || !r2->abs_error) return 1;

    // The references, independent of Numerixx: E'(1) and E''(1) for e = 0.5, from mpmath at 60 digits (E solved there
    // to 60 digits: E = 1.49870113351784831405798549726; Codex's 85-digit solve agrees). Each literal is within half an
    // ulp once rounded to a double: 1.1e-16 for E', 5.6e-17 for E''.
    double const d1 = 1.0373620218936458704990245805;
    double const d2 = -0.556713032668587779463634931485;
    // The page's closed forms, at the E the root solve gives, agree with them: that E is within about 1.6e-15 of the
    // root (the roots tutorial's bound), which moves 1 / (1 - e cos E) by |E''| 1.6e-15 = 9e-16 and -e sin E E'^3 by
    // |E'''| 1.6e-15, below 1.4e-15 (|E'''| = 0.855), plus a few ulps of rounding in cos, sin and the products: 1e-14
    // covers both with room.
    if (!(std::abs(exact - d1) < 1e-14 && std::abs(exact2 - d2) < 1e-14)) return 9;
    // "The first ten [decimals] are 1.0373620219": the printed slope and E' both round to it.
    if (!(std::format("{:.10f}", slope->value) == "1.0373620219" && std::format("{:.10f}", d1) == "1.0373620219"))
        return 7;
    // The central difference: "below 1e-9" relative. eccentric's values are E to within about 1.6e-15 (above), so the
    // rounding part is at most 2 * 1.6e-15 / 2h with h = eps^(1/3) = 6.06e-6: 2.64e-10 absolute; the truncation part is
    // h^2 / 6 |E'''|, below 1.3e-11; plus 1.1e-16 for d1's rounding.
    if (!(std::abs(slope->value - d1) <= 2.64e-10 + 1.3e-11 + 1.1e-16 && std::abs(slope->value - d1) / d1 < 1e-9))
        return 2;
    // "Its default step is about 6e-6 here": eps^(1/3) max(|M|, 1) with M = 1, 6.06e-6.
    if (!(dv::default_step<dv::Central<2>>(M) > 6e-6 && dv::default_step<dv::Central<2>>(M) < 6.1e-6)) return 8;
    // Ridders: the estimate covers the actual error (allowing d1's rounding), and is below 1e-11.
    if (!(std::abs(r->value - d1) <= *r->abs_error + 1.1e-16 && *r->abs_error < 1e-11)) return 3;
    // "More than ten times smaller than the central difference's, for about a dozen evaluations" (measured: 12; the
    // range allows for the platforms' std::sin).
    if (!(std::abs(r->value - d1) < std::abs(slope->value - d1) / 10)) return 5;
    int calls = 0;
    (void)dv::derivative<dv::Ridders>(
        [&](double m) {
            ++calls;
            return eccentric(m);
        },
        M);
    if (!(calls >= 8 && calls <= 20)) return 6;
    // The second derivative: the estimate covers the actual error (allowing d2's rounding), and is below 1e-10
    // (measured: 2.2e-11).
    if (!(std::abs(r2->value - d2) <= *r2->abs_error + 5.6e-17 && *r2->abs_error < 1e-10)) return 4;
    return 0;
}
