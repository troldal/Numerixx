// Compile-fail (DESIGN §3.3, §6.8): a bare number as a bracketing solver's stop criterion. It was a class template
// argument deduction failure without a reason; a deleted constructor now says what to write.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::bisection { nxx::width_tol { 1e-10 } };
#else
    auto s = r::bisection { 1e-10 };
#endif
    (void)s;
}
