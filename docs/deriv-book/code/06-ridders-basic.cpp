// E9 (Chapter 6): Ridders' extrapolation on the thread, cp = dH/dT at T0, with a wrapper that counts the
// evaluations of fn. Prints the value, the error estimate and the count; the actual error is at most 6 times the
// estimate (the tests' criterion; measured ratio 0.74) and the count at most 64 (the level cap for f'; measured 14).
// Quoted with \snippet{06-ridders-basic}{ridders-cp}.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using namespace nxx::deriv;
    double const T0 = book::thread_T0;
    // clang-format off
// <ridders-cp>
// Ridders from the default initial step T0 / 32, with fn
// wrapped to count how often it is called.
int calls = 0;
auto const H = [&calls](double T) {
    ++calls;
    return book::shomate_H(book::nitrogen, T);
};
auto const r = derivative<Ridders>(H, T0);
if (r)
    std::println("cp = {:.12f} +- {:.1e}, {} evaluations",
        r->value, *r->abs_error, calls);
// </ridders-cp>
    // clang-format on
    if (!r) return 1;
    double const cp = book::shomate_cp(book::nitrogen, T0);
    double const actual = std::abs(r->value - cp);
    std::println("exact {:.12f}, actual error {:.1e}, ratio {:.2f}", cp, actual, actual / *r->abs_error);
    bool const ok = actual <= 6 * *r->abs_error && calls <= 64;
    return ok ? 0 : 1;
}
