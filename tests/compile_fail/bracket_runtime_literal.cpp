// Compile-fail (DESIGN §6.2, §9.1): the literal constructor of nxx::bracket with a run-time endpoint. It is consteval, so
// a run-time value is not a constant expression; run-time values go through bracket<T>::make or a braced {lo, hi}.
#include <numerixx/core.hpp>

int main(int argc, char**)
{
    double lo = static_cast<double>(argc) - 1.0;
#ifdef NUMERIXX_CF_CONTROL
    auto b = nxx::bracket<double>::make(lo, 2.0);
#else
    auto b = nxx::bracket { lo, 2.0 };
#endif
    (void)b;
}
