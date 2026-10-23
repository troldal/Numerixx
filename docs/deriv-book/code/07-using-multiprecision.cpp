// E14 (Chapter 7, section 7.6; values DS18): derivative() on a user-defined real type, Boost.Multiprecision's
// cpp_bin_float_50 (50 decimal digits, epsilon 2^-167, about 5e-51). exp(x) at x = 1, whose derivative is e: the
// sixth-order central rule with its default step, and Ridders with its error estimate. Built only when the target
// boost_multiprecision exists (PLAN.md, decision A2).

#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <print>

using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Ridders;

int main()
{
    // clang-format off
// <multiprecision>
using mp50 = boost::multiprecision::cpp_bin_float_50;
// -> mp50: make the result a number, not an expression
// template that refers to x (Boost's types may return one).
auto const f = [](mp50 const &x) -> mp50 { return exp(x); };
mp50 const x = 1;
auto const d = derivative<Central<6>>(f, x);
// A double step converts to mp50: h0 = 1/32.
auto const r = derivative<Ridders>(f, x, 0.03125);
if (!d || !r) return 1;
// Boost's operators do not see the conversion of an
// Estimate: write mp50(*r) or r->value.
mp50 const e = exp(mp50(1));
mp50 const error = abs(mp50(*r) - e);
std::println("Central<6>: {}", d->value.str(50));
std::println("Ridders:    {}", r->value.str(50));
std::println("error {} estimate {}", error.str(2),
    r->abs_error->str(2));
// </multiprecision>
    // clang-format on
    mp50 const rel_c6 = abs(mp50(d->value - e)) / e;
    std::println("Central<6> relative error {}", rel_c6.str(2));
    // Central<6> at its default step eps^(1/7): the model gives about eps^(6/7) = 4e-44 relative; the generator
    // measures 7.2e-45 (DS18). The text claims at most 1e-40.
    bool ok = rel_c6 <= mp50(1e-40);
    // Ridders: the actual error at most 6 times its estimate, the tests' criterion (DS18: 9.9e-48 against
    // 4.4e-47 relative).
    ok = ok && r->abs_error.has_value() && error <= 6 * *r->abs_error;
    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}
