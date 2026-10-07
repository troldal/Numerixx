// Compile-fail (DESIGN §3.3 tier A, §6.2): the absolute part alone as an x_tol. x_tol{abs_tolerance{a}} would read as
// a purely absolute test, but a part alone is not a criterion; x_tol has its own deletion, whose text names x_tol and
// both run-time paths. The only case that reaches the abs_tolerance path of the part-alone deletion and its guide
// (§12.21 item 5); rel_tolerance_as_tolerance and width_tol_relative_alone reach the rel_tolerance path. The variable
// is not constexpr, so that cl reports C2280 at the deleted declaration rather than only C2131 (§9).
#include <numerixx/core.hpp>

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto c = nxx::x_tol { nxx::abs_tolerance { 1e-10 }, nxx::rel_tolerance { 1e-8 } };
#else
    auto c = nxx::x_tol { nxx::abs_tolerance { 1e-10 } };
#endif
    (void)c;
}
