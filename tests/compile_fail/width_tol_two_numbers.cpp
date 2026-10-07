// Compile-fail (DESIGN §3.3 tier A, §6.2): a width_tol literal from two bare numbers. Which one is relative? The
// roles could be swapped silently, so the relative part is always named with nxx::rel_tolerance, and the deleted
// constructor says so (reached through its deletion guide, not a CTAD error). The variable is not constexpr: from a
// constexpr initializer cl reports only C2131, not C2280 at the deleted declaration (§9).
#include <numerixx/core.hpp>

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto c = nxx::width_tol { 1e-10, nxx::rel_tolerance { 1e-8 } };
#else
    auto c = nxx::width_tol { 1e-10, 1e-8 };
#endif
    (void)c;
}
