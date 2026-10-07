// Compile-fail (DESIGN §3.3 tier A, §6.2): a rel_tolerance where a tolerance is expected. The roles are distinct types,
// not interchangeable; a part alone reaches the reasoned deletion of x_tol.
#include <numerixx/core.hpp>

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    nxx::x_tol<double> c { nxx::tolerance<double> { 0.1 } };
#else
    nxx::x_tol<double> c { nxx::rel_tolerance<double> { 0.1 } };
#endif
    return c.abs() > 0.0 ? 0 : 1;
}
