// Compile-fail (DESIGN §6.6): a braced {lo, hi} bracket with a function that cannot take its scalar type gets the
// facade's reason, like a bracket<T> or a std::pair does.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    (void)r::brent {}([](double x) { return x * x - 2.0; }, { 1.0, 2.0 });
#else
    (void)r::brent {}([](const char*) { return 0.0; }, { 1.0, 2.0 });
#endif
}
