// Compile-fail (DESIGN §6.10, §9.1): first_of over solvers that were not curried with .on(input), called with f alone.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    (void)nxx::first_of(r::brent {}.on(nxx::bracket { 0.0, 2.0 }), r::bisection {}.on(nxx::bracket { 0.0, 2.0 }))(f);
#else
    (void)nxx::first_of(r::brent {}, r::bisection {})(f);
#endif
}
