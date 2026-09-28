// P2564 probe (DESIGN §6.2 FLAG), not a compile-fail case: legal C++23 that GCC, Clang and clang-cl compile and cl
// rejects with C7595. make_solver forwards a raw scalar into a literal-checked tolerance, so under P2564 (consteval
// escalation) its specialisation becomes an immediate function, callable with a constant argument. MSVC lacks P2564,
// which is why every Numerixx function that forwards a tolerance takes the refined type (nxx::tolerance<T>) instead.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

template<class T>
constexpr auto make_solver(T t)
{ return r::bisection { nxx::width_tol { t } }; }

int main()
{
    auto           f   = [](double x) { return x * x - 2.0; };
    constexpr auto s   = make_solver(1e-6);
    auto           res = s(f, nxx::bracket { 0.0, 2.0 });
    return res ? 0 : 1;
}
