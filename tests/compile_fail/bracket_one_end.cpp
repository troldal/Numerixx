// Compile-fail (DESIGN §6.6): a braced list with one element is not a bracket. It once bound to the two-element array
// overload as {x, 0}, and the solve succeeded on [0, x], a bracket nobody wrote.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    const auto f = [](double x) { return x * x - 0.25; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::brent {}(f, { 0.0, 1.0 });
#else
    (void)r::brent {}(f, { 1.0 });
#endif
}
