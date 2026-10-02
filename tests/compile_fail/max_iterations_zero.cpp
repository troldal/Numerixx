// Compile-fail (DESIGN §3.3 tier A, §6.2): a zero iteration budget literal.
#include <numerixx/core.hpp>

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    nxx::max_iterations m = 1;
#else
    nxx::max_iterations m = 0;
#endif
    return static_cast<int>(m.value());
}
