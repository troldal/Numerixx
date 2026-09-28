// Compile-fail (DESIGN §3.3 tier A, §6.2): a negative tolerance literal, through alias-template CTAD.
#include <numerixx/core.hpp>

#ifdef NUMERIXX_CF_CONTROL
constexpr nxx::tolerance t { 1e-8 };
#else
constexpr nxx::tolerance t { -1e-8 };
#endif

int main() { return t.value() > 0.0 ? 0 : 1; }
