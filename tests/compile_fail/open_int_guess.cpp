// Compile-fail (DESIGN §3.3 tier A, §9.1): an open method given an int guess. The deleted overload of the open facade
// says how to write it.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::secant {}(f, 1.0);
#else
    (void)r::secant {}(f, 1);
#endif
}
