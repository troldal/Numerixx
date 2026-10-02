// Compile-fail (DESIGN §6.6, §6.13): a braced list with one element is not a bracket for the solve facade either. With a
// fixed [2] it bound as {x, 0}, and solve(f, {1.0}) succeeded on [0, 1], a bracket nobody wrote.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    const auto f = [](double x) { return x * x - 0.25; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::solve(f, { 0.0, 1.0 });
#else
    (void)r::solve(f, { 1.0 });
#endif
}
