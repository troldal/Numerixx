// Compile-fail (DESIGN §3.3, §6.8): a bare number as an open method's stop criterion. It was a class template argument
// deduction failure without a reason; a deleted constructor now says what to write.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::secant { nxx::x_tol { 1e-10 } };
#else
    auto s = r::secant { 1e-10 };
#endif
    (void)s;
}
