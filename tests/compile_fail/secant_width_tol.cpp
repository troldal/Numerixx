// Compile-fail (DESIGN §3.3, §6.8): width_tol on an open method (secant), whose views have no enclosure.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::secant { nxx::x_tol { 1e-6 } };
#else
    auto s = r::secant { nxx::width_tol { 1e-6 } };
#endif
    (void)s;
}
