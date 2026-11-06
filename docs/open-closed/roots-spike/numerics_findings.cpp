// Reproductions of the numerics review's findings 1 and 2 on draft 1 (PHASE7_API.md, Decisions from the reviews),
// run against today's header and the spike's. Build with -DNDEBUG to see the results past today's assert.
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <format>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace rt = nxx::roots;
using bin50 = boost::multiprecision::cpp_bin_float_50;
using dec50 = boost::multiprecision::cpp_dec_float_50;

template<typename T>
std::string show(std::expected<rt::Root<T>, rt::error> const &r, T const &unit)
{
    if (!r) return std::format("error {}", r.error());
    std::ostringstream o;
    o << std::setprecision(6) << T(r->value / unit) << " units, abs_error " << T(r->abs_error / unit) << " units";
    return o.str();
}

// Finding 1: the bracket [min, 2.5 min] around 1.25 min, default criteria, in a type without subnormals.
template<typename T>
void finding1(char const *type)
{
    T const min = (std::numeric_limits<T>::min)();
    auto const fn = [min](T const &x) { return T(T(x / min) - T(1.25)); };
    rt::Bracket<T> const b{ min, T(T(2.5) * min) };
    std::cout << type << ", [min, 2.5 min], root 1.25 min; 2 tol(min) = 2 min:\n";
    std::cout << "  Bisection      " << show(rt::solve<rt::Bisection>(fn, b), min) << '\n';
    std::cout << "  Brent          " << show(rt::solve<rt::Brent>(fn, b), min) << '\n';
#if defined(NDEBUG) || defined(ALL_METHODS) // ALL_METHODS: the spike, whose fix removes the assert
    std::cout << "  Toms748        " << show(rt::solve<rt::Toms748>(fn, b), min) << '\n';
    std::cout << "  AndersonBjorck " << show(rt::solve<rt::AndersonBjorck>(fn, b), min) << '\n';
    std::cout << "  Ridders        " << show(rt::solve<rt::Ridders>(fn, b), min) << '\n';
#endif
}

int main()
{
    finding1<bin50>("cpp_bin_float_50");
    finding1<dec50>("cpp_dec_float_50");

    // Finding 2: Bisection's stop test where the midpoint is exact.
    auto const shifted = [](double x) { return x - 2.4; };
    std::cout << "x - 2.4 on [1, 3], rel_tol 0.5, Bisection: "
              << show(rt::solve<rt::Bisection>(shifted, rt::Bracket{ 1.0, 3.0 }, { .rel_tol = 0.5 }), 1.0) << '\n';
    auto const pole = [](double x) { return 1.0 / x; };
    std::cout << "1/x on [-1, 1], abs_tol 1, Bisection: "
              << show(rt::solve<rt::Bisection>(pole, rt::Bracket{ -1.0, 1.0 }, { .abs_tol = 1.0 }), 1.0) << '\n';
}
