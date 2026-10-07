// Compile-fail (DESIGN §3.3 tier A, §6.2): a mixed width_tol literal with a negative absolute part. The absolute part is
// an abs_tolerance, so its own literal check rejects it (finite and >= 0); width_tol{0.0, nxx::rel_tolerance{1e-8}}
// shows that 0 itself is legal.
#include <numerixx/core.hpp>

#ifdef NUMERIXX_CF_CONTROL
constexpr auto c = nxx::width_tol { 0.0, nxx::rel_tolerance { 1e-8 } };
#else
constexpr auto c = nxx::width_tol { -1e-10, nxx::rel_tolerance { 1e-8 } };
#endif

int main() { return c.rel() > 0.0 ? 0 : 1; }
