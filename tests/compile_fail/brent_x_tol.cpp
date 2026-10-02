// Compile-fail (DESIGN §3.3, §6.8, §10.2 criterion 7): x_tol as Brent's tolerance, which must be a width criterion.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::brent { nxx::width_tol { 1e-6 } };
#else
    auto s = r::brent { nxx::x_tol { 1e-6 } };
#endif
    (void)s;
}
