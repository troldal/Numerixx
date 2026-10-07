// Compile-fail (DESIGN §3.3 tier A, §6.2): a part alone as a width criterion. width_tol{rel} would read as a purely
// relative test, but one argument is always absolute; the deletion (keyed on is_tolerance_part_v, so
// width_tol{nxx::abs_tolerance{a}} gets it too) names both spellings and their run-time paths. The variable is not
// constexpr, so that cl reports C2280 at the deleted declaration rather than only C2131 (§9).
#include <numerixx/core.hpp>

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto c = nxx::width_tol { 0.0, nxx::rel_tolerance { 1e-8 } };
#else
    auto c = nxx::width_tol { nxx::rel_tolerance { 1e-8 } };
#endif
    (void)c;
}
