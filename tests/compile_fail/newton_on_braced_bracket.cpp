// Compile-fail (DESIGN §3.3, §6.6): .on({lo, hi}) on an open method (Newton with a derivative). As with the call
// operator, the braced list had no overload; a deleted .on overload now says what to write.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto df = [](double x) { return 2.0 * x; };
#ifdef NUMERIXX_CF_CONTROL
    (void)r::newton {}.with_derivative(df).on(1.0);
#else
    (void)r::newton {}.with_derivative(df).on({ 1.0, 2.0 });
#endif
}
