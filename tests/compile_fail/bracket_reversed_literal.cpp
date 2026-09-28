// Compile-fail (DESIGN §3.3 tier A, §6.2): an invalid bracket literal (lo > hi).
#include <numerixx/core.hpp>

#ifdef NUMERIXX_CF_CONTROL
constexpr auto b = nxx::bracket { 1.0, 2.0 };
#else
constexpr auto b = nxx::bracket { 2.0, 1.0 };
#endif

int main() { return b.lo() < b.hi() ? 0 : 1; }
