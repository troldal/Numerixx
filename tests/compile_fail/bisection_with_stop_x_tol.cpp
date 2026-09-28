// Compile-fail (DESIGN §3.3, §6.8, §10.2 criterion 7): x_tol on a bracketing solver, through the with_stop builder.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    (void)r::bisection {}.with_stop(nxx::width_tol { 1e-6 });
#else
    (void)r::bisection {}.with_stop(nxx::x_tol { 1e-6 });
#endif
}
