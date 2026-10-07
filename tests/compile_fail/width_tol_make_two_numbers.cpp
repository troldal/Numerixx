// Compile-fail (DESIGN §3.3 tier A, §6.2): width_tol<T>::make(a, b) from two run-time numbers. make(T, T) is deleted
// so that the roles cannot be swapped: the relative part goes through rel_tolerance<T>::make first.
#include <numerixx/core.hpp>

int main()
{
    volatile double a = 1e-10;
    volatile double b = 1e-8;
#ifdef NUMERIXX_CF_CONTROL
    const auto w = nxx::rel_tolerance<double>::make(b).and_then([&](auto rp) { return nxx::width_tol<double>::make(a, rp); });
#else
    const auto w = nxx::width_tol<double>::make(a, b);
#endif
    return w ? 0 : 1;
}
