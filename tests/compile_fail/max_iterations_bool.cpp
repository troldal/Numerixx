// Compile-fail (DESIGN §3.3 tier A, §6.2): a bool as an iteration budget.
#include <numerixx/core.hpp>

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    nxx::max_iterations m = 10;
#else
    nxx::max_iterations m = true;
#endif
    return static_cast<int>(m.value());
}
