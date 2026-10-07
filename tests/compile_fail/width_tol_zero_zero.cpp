// Compile-fail (DESIGN §3.3 tier A, §6.2): a mixed width_tol literal with abs == 0 and rel == 0, which never converges.
// The parts are valid on their own (abs_tolerance accepts 0, rel_tolerance accepts 0), so it is the joint invariant of
// the consteval constructor that rejects it. width_tol{0.0, nxx::rel_tolerance{1e-8}} (purely relative) is legal.
#include <numerixx/core.hpp>

#ifdef NUMERIXX_CF_CONTROL
constexpr auto c = nxx::width_tol { 0.0, nxx::rel_tolerance { 1e-8 } };
#else
constexpr auto c = nxx::width_tol { 0.0, nxx::rel_tolerance { 0.0 } };
#endif

int main() { return c.rel() > 0.0 ? 0 : 1; }
