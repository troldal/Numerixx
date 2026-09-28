// Compile-fail (DESIGN §3.3 tier A, §9.1): a bracketing solver given a guess. The deleted overload of the bracketing
// facade names the inputs it takes.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::bisection {}(f, nxx::bracket { 0.0, 2.0 });
#else
    (void)r::bisection {}(f, 1.0);
#endif
}
