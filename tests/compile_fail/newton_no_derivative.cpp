// Compile-fail (DESIGN §3.3, §7.2, §9.1): Newton called without a derivative source, on a callable without
// .derivative().
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    auto df = [](double x) { return 2.0 * x; };
    (void)r::newton {}.with_derivative(df)(f, 1.0);
#else
    (void)r::newton {}(f, 1.0);
#endif
}
