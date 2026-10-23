// derivative() calls fn exactly as its constraint std::invocable<F &, T const &> promises: through std::invoke, with a
// T const lvalue, for every stencil. Regression tests: the evaluator for more than two nodes used to pass a temporary,
// and every branch called fn(...) directly, which rejects member-function pointers.

#include <numerixx/deriv.hpp>

#include <doctest/doctest.h>

#include <compare>
#include <concepts>
#include <limits>
#include <type_traits>

using nxx::deriv::Backward;
using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Forward;

namespace {

    // Accepts lvalues only: a temporary argument would select the deleted overload and fail to compile.
    struct LvalueOnly
    {
        constexpr double operator()(double const &x) const { return x * x; }
        double operator()(double &&) const = delete;
    };

    // A minimal user-defined real type with a member function, so that &Real::square is a member-function pointer
    // that std::invocable<F &, Real const &> accepts.
    struct Real
    {
        double v = 0;

        constexpr Real() = default;
        constexpr Real(int i) : v(i) {}
        constexpr explicit Real(double d) : v(d) {}

        [[nodiscard]] constexpr Real square() const { return Real(v * v); }

        friend constexpr Real operator+(Real a, Real b) { return Real(a.v + b.v); }
        friend constexpr Real operator-(Real a, Real b) { return Real(a.v - b.v); }
        friend constexpr Real operator*(Real a, Real b) { return Real(a.v * b.v); }
        friend constexpr Real operator/(Real a, Real b) { return Real(a.v / b.v); }
        friend constexpr Real operator-(Real a) { return Real(-a.v); }
        friend constexpr auto operator<=>(Real, Real) = default;
    };

} // namespace

template<>
struct std::numeric_limits<Real> : std::numeric_limits<double>
{
    static constexpr Real epsilon() { return Real(std::numeric_limits<double>::epsilon()); }
    static constexpr Real max() { return Real(std::numeric_limits<double>::max()); }
};

namespace {

    // A real type without a default constructor: nxx::real does not require one, so Ridders' tableau must not either.
    struct NoDefault
    {
        double v;

        constexpr NoDefault(int i) : v(i) {}
        constexpr explicit NoDefault(double d) : v(d) {}

        friend constexpr NoDefault operator+(NoDefault a, NoDefault b) { return NoDefault(a.v + b.v); }
        friend constexpr NoDefault operator-(NoDefault a, NoDefault b) { return NoDefault(a.v - b.v); }
        friend constexpr NoDefault operator*(NoDefault a, NoDefault b) { return NoDefault(a.v * b.v); }
        friend constexpr NoDefault operator/(NoDefault a, NoDefault b) { return NoDefault(a.v / b.v); }
        friend constexpr NoDefault operator-(NoDefault a) { return NoDefault(-a.v); }
        friend constexpr auto operator<=>(NoDefault, NoDefault) = default;
    };

} // namespace

template<>
struct std::numeric_limits<NoDefault> : std::numeric_limits<double>
{
    static constexpr NoDefault epsilon() { return NoDefault(std::numeric_limits<double>::epsilon()); }
    static constexpr NoDefault max() { return NoDefault(std::numeric_limits<double>::max()); }
};

static_assert(nxx::real<NoDefault> && !std::is_default_constructible_v<NoDefault>);

// Ridders on the user types: x^2 at 3 with h = 1/2 is exact at level 1 (A[1][1] equals its parent), so 6 exactly.
static_assert(derivative<nxx::deriv::Ridders>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
static_assert(
    derivative<nxx::deriv::Ridders>([](NoDefault const &x) { return x * x; }, NoDefault(3), NoDefault(0.5))->value.v
    == 6.0);

// Higher orders on the user types: x^2 has f'' = 2, and the second difference is exact for it whatever the step (2 at
// every Ridders level, so the tableau stops on exact agreement); f''' and f'''' of x^2 are 0.
static_assert(derivative<Central<2>, 2>(&Real::square, Real(3), Real(0.5))->value.v == 2.0);
static_assert(derivative<nxx::deriv::Ridders, 2>(&Real::square, Real(3), Real(0.5))->value.v == 2.0);
static_assert(derivative<Central<4>, 3>(&Real::square, Real(3), Real(0.5))->value.v == 0.0);
static_assert(
    derivative<nxx::deriv::Ridders, 4>([](NoDefault const &x) { return x * x; }, NoDefault(3), NoDefault(0.5))->value.v
    == 0.0);

static_assert(std::invocable<LvalueOnly &, double const &> && !std::invocable<LvalueOnly &, double>);
static_assert(nxx::real<Real> && std::invocable<decltype(&Real::square) &, Real const &>);

// x^2 at x = 1 with h = 1/2: the rules of order >= 2 are exact (2), Backward<1> has the error -h f''/2 = -1/2.
static_assert(derivative<Central<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
static_assert(derivative<Central<4>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
static_assert(derivative<Forward<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
static_assert(derivative<Backward<1>>(LvalueOnly{}, 1.0, 0.5) == 1.5);

// x^2 at x = 3 with h = 1/2, through the member-function pointer: exact, 6.
static_assert(derivative(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
static_assert(derivative<Central<4>>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
static_assert(derivative<Forward<2>>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);

TEST_CASE("deriv: callables are invoked with a T const lvalue through std::invoke")
{
    // The compile-time cases above, at run time.
    CHECK(derivative<Central<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<Central<4>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<Forward<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<Backward<1>>(LvalueOnly{}, 1.0, 0.5) == 1.5);
    CHECK(derivative(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
    CHECK(derivative<Central<4>>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
    CHECK(derivative<Forward<2>>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);

    // Ridders on the user types, as in the static_asserts above.
    CHECK(derivative<nxx::deriv::Ridders>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
    CHECK(
        derivative<nxx::deriv::Ridders>([](NoDefault const &x) { return x * x; }, NoDefault(3), NoDefault(0.5))->value.v
        == 6.0);

    // The higher orders above.
    CHECK(derivative<Central<2>, 2>(&Real::square, Real(3), Real(0.5))->value.v == 2.0);
    CHECK(derivative<nxx::deriv::Ridders, 2>(&Real::square, Real(3), Real(0.5))->value.v == 2.0);
    CHECK(derivative<Central<4>, 3>(&Real::square, Real(3), Real(0.5))->value.v == 0.0);
    CHECK(derivative<nxx::deriv::Ridders, 4>([](NoDefault const &x) { return x * x; }, NoDefault(3), NoDefault(0.5))
              ->value.v
          == 0.0);

    // The default step on the user type: x^2 is exact for the central difference, so only rounding remains. At x = 3,
    // h = 3 eps^(1/3); fn's values give S eps |f| / h = 3 eps^(2/3) (S = 1, f = 9) and the points
    // S/2 eps |x f'| / h = 3 eps^(2/3) (x f' = 18), so about 6 eps^(2/3) = 2.2e-10, or 4e-11 relative to 6.
    auto const d = derivative(&Real::square, Real(3));
    REQUIRE(d.has_value());
    CHECK(d->value.v == doctest::Approx(6.0).epsilon(2e-10));
}
