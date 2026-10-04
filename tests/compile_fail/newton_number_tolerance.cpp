// Compile-fail (DESIGN §3.3, §6.8): a bare number as Newton's stop criterion. newton{1e-10} failed class template
// argument deduction with no reason; a deleted constructor now says what to write.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::newton { nxx::x_tol { 1e-10 } };
#else
    auto s = r::newton { 1e-10 };
#endif
    (void)s;
}
