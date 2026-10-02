// Compile-fail (DESIGN §6.6, §6.13): the solve facade needs a bracket, as the bracketing solvers do. A scalar guess
// used to fail with a bare "no matching function", where brent{}(f, 1.0) states the reason.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    const auto f = [](double x) { return x * x - 0.25; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::solve(f, { 0.0, 1.0 });
#else
    (void)r::solve(f, 1.0);
#endif
}
