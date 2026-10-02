// Compile-fail (DESIGN §3.3, §7.2): a derivative on the secant method, which is derivative-free by design (only newton
// uses a derivative source), so the builder would configure something that means nothing.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    constexpr auto df = [](double x) { return 2.0 * x; };
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::newton {}.with_derivative(df);
#else
    auto s = r::secant {}.with_derivative(df);
#endif
    (void)s;
}
