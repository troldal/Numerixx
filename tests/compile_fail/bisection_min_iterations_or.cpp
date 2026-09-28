// Compile-fail (DESIGN §3.4, §6.8): under ||, min_iterations can still stop a solver on its own, so the combination is
// a guard too; here through the constructor.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    (void)r::bisection { nxx::width_tol { 1e-6 } && nxx::min_iterations { 3 } };
#else
    (void)r::bisection { nxx::width_tol { 1e-6 } || nxx::min_iterations { 3 } };
#endif
}
