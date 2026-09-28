// Compile-fail (DESIGN §3.4, §6.8): min_iterations alone would report success after n iterations without testing
// accuracy, so a solver takes it only under && with a convergence test; here through the with_stop builder.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    (void)r::secant {}.with_stop(nxx::x_tol { 1e-10 } && nxx::min_iterations { 3 });
#else
    (void)r::secant {}.with_stop(nxx::min_iterations { 3 });
#endif
}
