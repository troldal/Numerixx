// Compile-fail (DESIGN §3.3, §6.8): step_tol on a bracketing solver, through its constructor. Like x_tol, it compares
// successive iterates.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::bisection { nxx::floored_width { 40 } };
#else
    auto s = r::bisection { nxx::step_tol<3, 5> {} };
#endif
    (void)s;
}
