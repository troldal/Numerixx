// Compile-fail (DESIGN §6.3): nxx::best on a search result, which succeeds with a sign_bracket and fails with a
// root_estimate, reaches the deleted sibling, whose reason says to read *r and r.error().best separately.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    const auto res = r::brent {}(f, nxx::bracket { 1.0, 2.0 });
#else
    const auto res = r::expand {}(f, nxx::bracket { 2.0, 2.5 });
#endif
    (void)nxx::best(res);
}
