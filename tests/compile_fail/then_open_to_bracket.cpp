// Compile-fail (DESIGN §6.10, §9.1): then(open method, bracketing solver). A bracketing solver cannot start from a point
// estimate; the then contract says so when it is called.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f  = [](double x) { return x * x - 2.0; };
    auto df = [](double x) { return 2.0 * x; };
#ifdef NUMERIXX_CF_CONTROL
    (void)nxx::then(r::bisection {}.on(nxx::bracket { 0.0, 2.0 }), r::newton {}.with_derivative(df))(f);
#else
    (void)nxx::then(r::newton {}.with_derivative(df).on(1.0), r::bisection {})(f);
#endif
}
