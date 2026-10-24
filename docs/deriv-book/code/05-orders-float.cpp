// E8 (Chapter 5, quoted in Chapter 8): a fourth derivative in float. The call selects a deprecated overload of
// derivative(), whose warning says that the result is mostly rounding noise; this file builds with warnings as
// errors, so the warning must be silenced at the call (the pragma), or for the whole build with
// NUMERIXX_NO_PRECISION_WARNING. The function is exp(x) at x = 1 (every derivative of exp is exp), the same call
// the data generator makes for the \valFloatFourth...AtOne macros of DS6, so the printed numbers are the book's.
// The example only asserts that a value comes back. Quoted with \snippet{05-orders-float}{float-fourth}.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using namespace nxx::deriv;
    // clang-format off
// <float-fourth>
auto const f = [](float x) { return std::exp(x); };
// f'''' in float selects a deprecated overload: its warning
// says the result is mostly rounding noise. Silence it at
// the call, as here, or for the whole build by defining
// NUMERIXX_NO_PRECISION_WARNING.
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(push)
#pragma warning(disable : 4996)
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
auto const d = derivative<Central<2>, 4>(f, 1.0f);
auto const r = derivative<Ridders, 4>(f, 1.0f);
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif
// </float-fourth>
    // clang-format on
    double const exact = std::exp(1.0); // every fourth derivative of exp is exp
    if (!d || !r) return 1;
    auto const rel = [exact](float value) { return std::abs(static_cast<double>(value) - exact) / exact; };
    std::println("exp'''' (1) exact    {:.7f}", exact);
    std::println("Central<2>, float    {:.7f}  rel. error {:.1e}", d->value, rel(d->value));
    std::println("Ridders, float       {:.7f}  rel. error {:.1e}  estimate {:.1e} ({:.1e} relative)", r->value,
        rel(r->value), *r->abs_error, static_cast<double>(*r->abs_error) / exact);
    return 0;
}
