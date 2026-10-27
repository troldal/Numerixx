// nxx::deriv::derivative_of (docs/roots/PHASE3_API.md): a derivative as a callable. Its call operator forwards to the
// implementation behind derivative<M, N>(fn, x) without branching on M, N or T, so the tests compare the two bit for
// bit on a few methods and orders, and check what is its own: the copy of fn, its constness and the constraint.

#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>

using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::derivative_of;
using nxx::deriv::Ridders;

namespace {

    constexpr auto cube = [](double x) { return x * x * x - 2.0 * x; };

    constexpr double cube_free(double x) { return x * x * x - 2.0 * x; }

    // A function object whose parameter can be changed after derivative_of has copied it.
    struct Scaled
    {
        double k;
        constexpr double operator()(double x) const { return k * x * x; }
    };

    // Callable with float only: a double argument selects the deleted template.
    struct OnlyFloat
    {
        float operator()(float x) const { return x; }
        template<typename U>
        float operator()(U) const = delete;
    };

    auto outlived()
    {
        double const k = 3.0;
        return derivative_of([k](double x) { return k * x * x; }); // the lambda, and k with it, is copied
    }

} // namespace

// The same result as derivative<M, N>(fn, x), bit for bit, the whole std::expected: the order passes through (an
// integer or an order), and so does Ridders' error estimate.
static_assert(derivative_of(cube)(1.5) == derivative(cube, 1.5));
static_assert(derivative_of<Central<2>, nxx::deriv::order::second>(cube)(1.5)
              == derivative<Central<2>, nxx::deriv::order::second>(cube, 1.5));
static_assert(derivative_of<Ridders>(cube)(1.5) == derivative<Ridders>(cube, 1.5));
static_assert(derivative_of<Ridders>(cube)(1.5)->abs_error.has_value());
static_assert(derivative_of(cube_free)(1.5) == derivative(cube_free, 1.5)); // a free function, stored as a pointer
// A failure is the same error: x not finite.
static_assert(derivative_of(cube)(std::numeric_limits<double>::infinity())
              == std::unexpected(nxx::deriv::error::non_finite));

// The copy: a struct changed after derivative_of does not change the callable's result; std::ref shares it.
constexpr bool copies()
{
    Scaled s{ 2.0 };
    auto const copied = derivative_of(s);
    auto const shared = derivative_of(std::ref(s));
    s.k = 3.0;
    return copied(1.0) == derivative(Scaled{ 2.0 }, 1.0) && shared(1.0) == derivative(Scaled{ 3.0 }, 1.0);
}
static_assert(copies());

// A lambda capturing a local by reference, in constant evaluation: the reference is copied, the local shared.
constexpr bool counts()
{
    int n = 0;
    auto const d = derivative_of([&n](double x) {
        ++n;
        return x * x;
    });
    return d(1.0).has_value() && n == 2; // Central<2>'s two nodes
}
static_assert(counts());

// Constness as std::bind_front: a mutable fn's copy is callable through a non-const callable, not a const one.
using Mutable = decltype(derivative_of([k = 1.0](double x) mutable {
    k += 1.0;
    return k * x;
}));
static_assert(std::is_invocable_v<Mutable &, double const &>);
static_assert(!std::is_invocable_v<Mutable const &, double const &>);
static_assert(std::is_invocable_v<decltype(derivative_of(cube)) const &, double const &>);

// The constraint: a method and an order derivative accepts; a T fn accepts (the table's integers: test_callables).
template<typename M, auto N>
constexpr bool makes = requires { derivative_of<M, N>(cube); };
static_assert(makes<Central<2>, 1> && makes<Ridders, 4> && makes<Central<2>, nxx::deriv::order::third>);
static_assert(!makes<double, 1> && !makes<Central<2>, 0> && !makes<Central<2>, 5>);
static_assert(!std::is_invocable_v<decltype(derivative_of(cube)) const &, std::string const &>);
static_assert(!std::is_invocable_v<decltype(derivative_of(OnlyFloat{})) const &, double const &>);
static_assert(std::is_invocable_v<decltype(derivative_of(OnlyFloat{})) const &, float const &>);

TEST_CASE("deriv:of: the same result as derivative")
{
    CHECK(derivative_of(cube)(1.5) == derivative(cube, 1.5));
    CHECK(derivative_of<Central<2>, nxx::deriv::order::second>(cube)(1.5)
          == derivative<Central<2>, nxx::deriv::order::second>(cube, 1.5));
    CHECK(derivative_of<Ridders>(cube)(1.5) == derivative<Ridders>(cube, 1.5));

    // libm: the same calls in the same order, so the same values.
    auto const lg = [](double x) { return std::log(x); };
    CHECK(derivative_of<Ridders>(lg)(2.0) == derivative<Ridders>(lg, 2.0));
    // A failure: the stencil at 0 evaluates log(-h), NaN.
    CHECK(derivative_of(lg)(0.0) == std::unexpected(nxx::deriv::error::non_finite));
    CHECK(derivative(lg, 0.0) == std::unexpected(nxx::deriv::error::non_finite));

    using bin50 = boost::multiprecision::cpp_bin_float_50;
    auto const big = [](bin50 const &x) { return bin50(x * x * x); };
    CHECK(derivative_of(big)(bin50(1.5)) == derivative(big, bin50(1.5)));

    // A float second derivative through derivative_of: no deprecation warning (this file builds with warnings as
    // errors), and the same value as derivative's, which warns.
    auto const fc = [](float x) { return x * x * x; };
    auto const d2 = derivative_of<Central<2>, 2>(fc)(1.5f);
    REQUIRE(d2.has_value());
    CHECK(d2->value == doctest::Approx(9.0).epsilon(1e-2)); // float's f'' by Central<2>: about 1e-3 relative
}

TEST_CASE("deriv:of: the copy of fn")
{
    CHECK(copies());
    CHECK(counts());

    auto const d = outlived(); // the lambda's scope has ended
    auto const r = d(2.0);
    REQUIRE(r.has_value());
    CHECK(r->value == doctest::Approx(12.0).epsilon(1e-9)); // 6x at 2; Central<2> on a quadratic: rounding only

    // A mutable fn through a non-const callable: its copy's state changes, the original's does not.
    int n = 0;
    auto counter = [&n, calls = 0](double x) mutable {
        ++calls;
        n = calls;
        return x;
    };
    auto copy = derivative_of(counter);
    CHECK(copy(1.0).has_value());
    CHECK(copy(1.0).has_value());
    CHECK(n == 4); // the copy's own count, two nodes per call
    counter(1.0);
    CHECK(n == 1); // the original's count was never touched
}
