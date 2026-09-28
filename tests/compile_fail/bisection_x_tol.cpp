// Compile-fail (DESIGN §3.3, §6.8, §10.2 criterion 7): x_tol on a bracketing solver, through its constructor. x_tol
// compares successive iterates, which says nothing about an enclosure.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::bisection { nxx::width_tol { 1e-6 } };
#else
    auto s = r::bisection { nxx::x_tol { 1e-6 } };
#endif
    (void)s;
}
