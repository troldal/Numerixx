// Compile-fail (DESIGN §6.6, §7.2): a validated tolerance is not a criterion. bisection has no guide for it: the
// implicit guide from its defaulted Opt reaches the widened bare-number deletion, as it does for a bare number.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    if (auto tol = nxx::tolerance<double>::make(1e-10)) {
#ifdef NUMERIXX_CF_CONTROL
        auto s = r::bisection { nxx::width_tol { *tol } };
#else
        auto s = r::bisection { *tol };
#endif
        (void)s;
    }
}
