// Compile-fail (DESIGN §6.6, §7.2): with_stop on a searcher, which has no configurable stop criterion.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    (void)r::expand {}.with_budget(10);
#else
    (void)r::expand {}.with_stop(nxx::never {});
#endif
}
