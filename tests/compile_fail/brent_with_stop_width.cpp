// Compile-fail (DESIGN §6.8): a width criterion given to brent through with_stop. brent's intrinsic test applies its
// own tolerance and runs before the stop criterion, so with_stop cannot tighten it: it reported stop_reason::criterion
// at brent's default tolerance (width 6.66e-16 for width_tol{1e-20} on x^2 - 2 over [1, 2]). The width criterion
// goes to the constructor.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    (void)r::brent { nxx::width_tol { 1e-12 } };
#else
    (void)r::brent {}.with_stop(nxx::width_tol { 1e-12 });
#endif
}
