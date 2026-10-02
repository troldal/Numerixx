// Compile-fail (DESIGN §3.3, §6.9, D20): a projection on a bracketing solver. Projection applies to open methods; a
// bracketing method keeps every iterate inside its bracket, so the builder would configure something that means nothing.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::secant {}.with_projection(r::clamp_to { 0.0, 2.0 });
#else
    auto s = r::bisection {}.with_projection(r::clamp_to { 0.0, 2.0 });
#endif
    (void)s;
}
