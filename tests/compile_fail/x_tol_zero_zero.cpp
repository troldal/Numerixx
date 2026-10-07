// Compile-fail (DESIGN §3.3 tier A, §6.2): a mixed x_tol literal with abs == 0 and rel == 0, which never converges.
// x_tol{0.0, nxx::rel_tolerance{1e-8}} (purely relative) is legal.
#include <numerixx/core.hpp>

#ifdef NUMERIXX_CF_CONTROL
constexpr auto c = nxx::x_tol { 0.0, nxx::rel_tolerance { 1e-8 } };
#else
constexpr auto c = nxx::x_tol { 0.0, nxx::rel_tolerance { 0.0 } };
#endif

int main() { return c.rel() > 0.0 ? 0 : 1; }
