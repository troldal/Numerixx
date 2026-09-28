// Self-test of the compile-fail harness (cmake/NumerixxCompileFail.cmake). The control build must compile; the
// negative build must fail, and on GCC and Clang its first error must carry the reason below. Phase 1 adds the
// real illegal-state cases (DESIGN §9.1).
#include <numerixx/core.hpp>

#ifdef NUMERIXX_CF_CONTROL
static_assert(nxx::version.major == NUMERIXX_VERSION_MAJOR);
#else
static_assert(nxx::version.major < 0, "compile-fail harness self-test");
#endif
