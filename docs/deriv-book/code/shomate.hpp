// The thread example of the book: the molar enthalpy of nitrogen gas in the Shomate form of the NIST Chemistry
// WebBook, with its closed-form derivatives up to order 6. The data generator (tools/generate_data.cpp) and the
// examples in code/ include this one file, so the coefficients exist in exactly one place.
//
// Source of the equations and the coefficients: P. J. Linstrom and W. G. Mallard (eds.), NIST Chemistry WebBook,
// NIST Standard Reference Database Number 69, National Institute of Standards and Technology, Gaithersburg MD,
// https://doi.org/10.18434/T4D303; nitrogen (CAS 7727-37-9), "Gas phase thermochemistry data", Shomate equation,
// https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1, retrieved 2026-10-23. The page states "Data last
// reviewed in March, 1977; New parameter fit January 2009" and names as the source of the data M. W. Chase, Jr.,
// NIST-JANAF Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data, Monograph 9, 1998 (nitrogen: Part II,
// p. 1621). The coefficients are the WebBook's fit to those tables; the tables themselves hold no coefficients.
// Every value below was verified against the page source by the book's evidence reviewer on 2026-10-23.
//
// Licence: the coefficients in nitrogen_low, nitrogen and nitrogen_high are NIST Standard Reference Data,
// (c) U.S. Secretary of Commerce, reproduced with attribution; they are not covered by the MIT licence of this
// repository, which applies to the code around them.
//
// The WebBook's equations, with t = T / 1000 and T in kelvin:
//     Cp                 = A + B t + C t^2 + D t^3 + E / t^2                      J/(mol K)
//     H - H(298.15 K)    = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F - H       kJ/mol
// This header returns H in J/mol (1000 times the WebBook's kJ/mol), so that dH/dT is Cp in J/(mol K) exactly: the
// factor 1000 of the unit cancels the factor 1/1000 of dt/dT. H is the enthalpy relative to 298.15 K; the
// coefficient named H (0.0 for nitrogen) is a constant of the fit, not the enthalpy.
//
// Units of the derivatives: d^n H / dT^n is in J/(mol K^n).

#pragma once

#include <numerixx/math.hpp>

#include <cassert>

namespace book {

// clang-format off
// Snippet lines stay within 60 characters (code/README.md).
// <shomate>
// One range of a Shomate fit: the WebBook's eight
// coefficients and the temperatures the fit is valid on.
struct Shomate
{
    double A, B, C, D, E, F, G, H;
    double T_min, T_max; // K
};

// Nitrogen (N2), NIST Chemistry WebBook, fit of January
// 2009 to Chase 1998; retrieved 2026-10-23. Three ranges
// that share their end points.
inline constexpr Shomate nitrogen_low{
    28.98641, 1.853978, -9.647459, 16.63537, 0.000117,
    -8.671914, 226.4168, 0.0, 100.0, 500.0 };
inline constexpr Shomate nitrogen{
    19.50583, 19.88705, -8.598535, 1.369784, 0.527601,
    -4.935202, 212.3900, 0.0, 500.0, 2000.0 };
inline constexpr Shomate nitrogen_high{
    35.51872, 1.128728, -0.196103, 0.014662, -4.553760,
    -18.97091, 224.9810, 0.0, 2000.0, 6000.0 };
// </shomate>

// <shomate-H>
// H(T) - H(298.15 K) in J/mol for T in K: the WebBook's
// kJ/mol expression times 1000. Valid on [T_min, T_max];
// evaluable anywhere above 0 K.
template<nxx::real Real>
constexpr Real shomate_H(Shomate const &c, Real const &T)
{
    Real const t = T / Real(1000);
    return Real(Real(1000)
        * (Real(c.A) * t + Real(c.B) * t * t / Real(2)
            + Real(c.C) * t * t * t / Real(3)
            + Real(c.D) * t * t * t * t / Real(4)
            - Real(c.E) / t + Real(c.F) - Real(c.H)));
}
// </shomate-H>

// <shomate-cp>
// cp(T) = dH/dT in J/(mol K): the WebBook's heat-capacity
// expression.
template<nxx::real Real>
constexpr Real shomate_cp(Shomate const &c, Real const &T)
{
    Real const t = T / Real(1000);
    return Real(Real(c.A) + Real(c.B) * t
        + Real(c.C) * t * t + Real(c.D) * t * t * t
        + Real(c.E) / (t * t));
}
// </shomate-cp>
// clang-format on

    // d^n H / dT^n for n = 0 to 6, in J/(mol K^n): the answer key for every derivative the book measures. Each
    // derivative with respect to T is the derivative with respect to t divided by 1000; from the fifth on only the
    // E/t term is left.
    template<nxx::real Real>
    constexpr Real shomate_derivative(Shomate const &c, Real const &T, int n)
    {
        assert(n >= 0 && n <= 6);
        Real const t = T / Real(1000);
        Real const k = Real(1) / Real(1000); // dt/dT
        Real const t2 = t * t;
        Real const t3 = t2 * t;
        Real const E = Real(c.E);
        switch (n) {
            case 0: return book::shomate_H(c, T);
            case 1: return book::shomate_cp(c, T);
            case 2: return Real(k * (Real(c.B) + Real(2) * Real(c.C) * t + Real(3) * Real(c.D) * t2 - Real(2) * E / t3));
            case 3: return Real(k * k * (Real(2) * Real(c.C) + Real(6) * Real(c.D) * t + Real(6) * E / (t3 * t)));
            case 4: return Real(k * k * k * (Real(6) * Real(c.D) - Real(24) * E / (t3 * t2)));
            case 5: return Real(k * k * k * k * (Real(120) * E / (t3 * t3)));
            default: return Real(k * k * k * k * k * (Real(-720) * E / (t3 * t3 * t)));
        }
    }

    // The thread's point: 1000 K, interior to the 500 to 2000 K range (its lower end, 500 K, is the domain edge of
    // Chapters 4 and 7).
    inline constexpr double thread_T0 = 1000.0;

    // The enthalpy of nitrogen across the junction of the two lower ranges at 500 K, each range used where it is
    // valid: the lower fit below 500 K, the thread's fit from 500 K on (the WebBook does not say which range applies
    // exactly at a shared end point; this header gives the end point to the upper range). The two fits do not meet
    // exactly: at 500 K, H jumps by about 0.77 J/mol and cp by about -7.5e-5 J/(mol K), and dcp/dT falls by 17 %
    // (the evidence reviewer's junction check, 2026-10-23). This is the range break of the smoothness pitfall
    // (Chapter 7; E19, DS20): the function is evaluable everywhere, but not smooth at 500 K.
    inline constexpr double junction_T = 500.0;

    template<nxx::real Real>
    constexpr Real nitrogen_H_joined(Real const &T)
    { return T < Real(junction_T) ? book::shomate_H(nitrogen_low, T) : book::shomate_H(nitrogen, T); }

    template<nxx::real Real>
    constexpr Real nitrogen_derivative_joined(Real const &T, int n)
    { return T < Real(junction_T) ? book::shomate_derivative(nitrogen_low, T, n) : book::shomate_derivative(nitrogen, T, n); }

} // namespace book
