// Compile-fail (DESIGN §6.6): a pointer is not a bracket, and .on says so (the catch-all takes a forwarding reference,
// so a C array still reaches the array overloads instead of decaying to a pointer).
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    const double ends[2] = { 1.0, 2.0 };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::bisection {}.on(ends);
#else
    const double* p = ends;
    (void)r::bisection {}.on(p);
#endif
}
