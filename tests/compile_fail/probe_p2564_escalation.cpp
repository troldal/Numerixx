// P2564 probe (DESIGN §6.2 FLAG), not a compile-fail case: legal C++23 that GCC, Clang and clang-cl compile and cl
// rejects with C7595. make_solver forwards a raw scalar into a literal-checked tolerance, so under P2564 (consteval
// escalation) its specialisation becomes an immediate function, callable with a constant argument. MSVC lacks P2564,
// which is why every Numerixx function that forwards a tolerance takes the refined type (nxx::tolerance<T>) instead.
// make_mixed forwards refined parts into the role-typed mixed constructor, which is consteval itself, so it escalates
// too and fails on cl the same way: generic code that forwards a mixed tolerance calls width_tol<T>::make instead.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

template<class T>
constexpr auto make_solver(T t)
{ return r::bisection { nxx::width_tol { t } }; }

template<class T>
constexpr auto make_mixed(nxx::abs_tolerance<T> a, nxx::rel_tolerance<T> rel)
{ return r::bisection { nxx::width_tol { a, rel } }; }

int main()
{
    auto           f   = [](double x) { return x * x - 2.0; };
    constexpr auto s   = make_solver(1e-6);
    constexpr auto m   = make_mixed(nxx::abs_tolerance { 1e-6 }, nxx::rel_tolerance { 1e-8 });
    auto           res = s(f, nxx::bracket { 0.0, 2.0 });
    auto           mix = m(f, nxx::bracket { 0.0, 2.0 });
    return res && mix ? 0 : 1;
}
