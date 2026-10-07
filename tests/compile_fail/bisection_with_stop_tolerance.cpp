// Compile-fail (DESIGN §6.6): with_stop given a validated tolerance. It reached the criterion catch-all and got the false
// reason "this criterion does not apply to this solver"; a deleted sibling for everything that is not a criterion now
// names the tests in x and what each one bounds.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    if (auto tol = nxx::tolerance<double>::make(1e-10)) {
#ifdef NUMERIXX_CF_CONTROL
        (void)r::bisection {}.with_stop(nxx::width_tol { *tol });
#else
        (void)r::bisection {}.with_stop(*tol);
#endif
    }
}
