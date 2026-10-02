// Compile-fail (DESIGN §6.10, §9.1): an any_solver built from a curried solver with a different result type (a searcher
// returns a sign_bracket, not a root_estimate).
#include <numerixx/core/any_solver.hpp>
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

namespace
{
    using fn_t     = double (*)(double);
    using solver_t = nxx::any_solver<fn_t, r::root_estimate<double>>;
}    // namespace

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    solver_t s = r::bisection {}.on(nxx::bracket { 0.0, 1.0 });
#else
    solver_t s = r::expand {}.on(nxx::bracket { 0.0, 1.0 });
#endif
    (void)s;
}
