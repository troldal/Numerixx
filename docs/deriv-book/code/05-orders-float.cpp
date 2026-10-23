// E8 (Chapter 5): a fourth derivative in float. The call selects a deprecated overload of derivative(), whose
// warning says that the result is mostly rounding noise; this file builds with warnings as errors, so the warning
// must be silenced at the call (the pragma), or for the whole build with NUMERIXX_NO_PRECISION_WARNING. The example
// only asserts that a value comes back; the text quotes DS6 for the accuracy. Quoted with
// \snippet{05-orders-float}{float-fourth}.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using namespace nxx::deriv;
    // clang-format off
// <float-fourth>
auto const f = [](float x) { return std::sin(x); };
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
auto const d = derivative<Central<2>, 4>(f, 0.5f);
auto const r = derivative<Ridders, 4>(f, 0.5f);
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif
// </float-fourth>
    // clang-format on
    float const exact = std::sin(0.5f); // every fourth derivative of sin is sin
    if (!d || !r) return 1;
    std::println("sin'''' (0.5) exact  {:.7f}", exact);
    std::println("Central<2>, float    {:.7f}  rel. error {:.1e}", d->value, std::abs(d->value - exact) / exact);
    std::println("Ridders, float       {:.7f}  rel. error {:.1e}  estimate {:.1e}", r->value,
        std::abs(r->value - exact) / exact, *r->abs_error);
    return 0;
}
