// Compile-fail (DESIGN §6.10, §9.1): first_of over alternatives with different result types (a bracketing solver and a
// searcher).
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    (void)nxx::first_of(r::bisection {}.on(nxx::bracket { 0.0, 2.0 }), r::brent {}.on(nxx::bracket { 0.0, 1.0 }))(f);
#else
    (void)nxx::first_of(r::bisection {}.on(nxx::bracket { 0.0, 2.0 }), r::expand {}.on(nxx::bracket { 0.0, 1.0 }))(f);
#endif
}
