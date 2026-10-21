// derivative() calls fn exactly as its constraint std::invocable<F &, T const &> promises: through std::invoke, with a
// T const lvalue, for every stencil. Regression tests: the evaluator for more than two nodes used to pass a temporary,
// and every branch called fn(...) directly, which rejects member-function pointers.

#include <numerixx/deriv.hpp>

#include <doctest/doctest.h>

#include <compare>
#include <concepts>
#include <limits>

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

static_assert(std::invocable<LvalueOnly &, double const &> && !std::invocable<LvalueOnly &, double>);
static_assert(nxx::real<Real> && std::invocable<decltype(&Real::square) &, Real const &>);

// x^2 at x = 1 with h = 1/2: the rules of order >= 2 are exact (2), Backward<1> has the error -h f''/2 = -1/2.
static_assert(derivative<Central<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
static_assert(derivative<Central<4>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
static_assert(derivative<Forward<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
static_assert(derivative<Backward<1>>(LvalueOnly{}, 1.0, 0.5) == 1.5);

// x^2 at x = 3 with h = 1/2, through the member-function pointer: exact, 6.
static_assert(derivative(&Real::square, Real(3), Real(0.5))->v == 6.0);
static_assert(derivative<Central<4>>(&Real::square, Real(3), Real(0.5))->v == 6.0);
static_assert(derivative<Forward<2>>(&Real::square, Real(3), Real(0.5))->v == 6.0);

TEST_CASE("deriv: callables are invoked with a T const lvalue through std::invoke")
{
    // The compile-time cases above, at run time.
    CHECK(derivative<Central<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<Central<4>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<Forward<2>>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<Backward<1>>(LvalueOnly{}, 1.0, 0.5) == 1.5);
    CHECK(derivative(&Real::square, Real(3), Real(0.5))->v == 6.0);
    CHECK(derivative<Central<4>>(&Real::square, Real(3), Real(0.5))->v == 6.0);
    CHECK(derivative<Forward<2>>(&Real::square, Real(3), Real(0.5))->v == 6.0);

    // The default step on the user type: x^2 is exact for the central difference, so only rounding remains. At x = 3,
    // h = 3 eps^(1/3); fn's values give S eps |f| / h = 3 eps^(2/3) (S = 1, f = 9) and the points
    // S/2 eps |x f'| / h = 3 eps^(2/3) (x f' = 18), so about 6 eps^(2/3) = 2.2e-10, or 4e-11 relative to 6.
    auto const d = derivative(&Real::square, Real(3));
    REQUIRE(d.has_value());
    CHECK(d->v == doctest::Approx(6.0).epsilon(2e-10));
}
