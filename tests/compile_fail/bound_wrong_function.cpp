// Compile-fail (DESIGN §6.6): a curried solver given a function it cannot take. bound::operator() was constrained with no
// reason (Clang: "no matching function for call to object of type 'bound<...>'"); its deleted sibling now says where the
// reason is.
#include <numerixx/roots.hpp>

#include <utility>

namespace r = nxx::roots;

int main()
{
    const auto curried = r::brent {}.on(std::pair { 1.0, 2.0 });
#ifdef NUMERIXX_CF_CONTROL
    (void)curried([](double x) { return x * x - 2.0; });
#else
    (void)curried([](const char*) { return 0.0; });
#endif
}
