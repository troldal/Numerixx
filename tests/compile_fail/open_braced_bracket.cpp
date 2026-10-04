// Compile-fail (DESIGN §3.3, §6.6): an open method given a braced {lo, hi}. open_facade had no overload for a braced
// list or a C array, so the call failed with a bare "no matching function"; a deleted overload now says what to write.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::secant {}(f, 1.0);
#else
    (void)r::secant {}(f, { 1.0, 2.0 });
#endif
}
