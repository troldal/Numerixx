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

namespace {

    // A real type whose numeric_limits is written member by member, with exactly the members nxx::real requires and
    // none inherited from numeric_limits<double> (U-5): every rule and Ridders must compile and run with it.
    struct Plain
    {
        double v;

        constexpr Plain(int i) : v(i) {}
        constexpr explicit Plain(double d) : v(d) {}

        friend constexpr Plain operator+(Plain a, Plain b) { return Plain(a.v + b.v); }
        friend constexpr Plain operator-(Plain a, Plain b) { return Plain(a.v - b.v); }
        friend constexpr Plain operator*(Plain a, Plain b) { return Plain(a.v * b.v); }
        friend constexpr Plain operator/(Plain a, Plain b) { return Plain(a.v / b.v); }
        friend constexpr Plain operator-(Plain a) { return Plain(-a.v); }
        friend constexpr auto operator<=>(Plain, Plain) = default;
    };

    constexpr Plain plain_square(Plain const &x) { return x * x; }

} // namespace

template<>
struct std::numeric_limits<Plain>
{
    static constexpr bool is_specialized = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_exact = false;
    static constexpr int radix = 2;
    static constexpr int digits = 53;
    static constexpr Plain epsilon() { return Plain(std::numeric_limits<double>::epsilon()); }
    static constexpr Plain max() { return Plain((std::numeric_limits<double>::max)()); }
    static constexpr Plain denorm_min() { return Plain(std::numeric_limits<double>::denorm_min()); }
};

static_assert(nxx::real<Plain>);

// x^2 at 3 with h = 1/2: exact for every rule of order >= 2 and for Ridders (6; the second derivative is 2, the
// third 0).
static_assert(derivative(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
static_assert(derivative<Forward<2>>(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
static_assert(derivative<nxx::deriv::Ridders>(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
static_assert(derivative<Central<2>, 2>(plain_square, Plain(3), Plain(0.5))->value.v == 2.0);
static_assert(derivative<Central<4>, 3>(plain_square, Plain(3), Plain(0.5))->value.v == 0.0);

namespace {

    // Types that declare few significant binary digits, numeric_limits written member by member (U-10): D = 11 is half
    // precision, D = 2 holds the integers 0 to 3. A template, so that its operators, which only the constraints
    // inspect in most cases, are not reported as unused.
    template<int D>
    struct Narrow
    {
        double v = 0;

        constexpr Narrow(int i) : v(i) {}
        constexpr explicit Narrow(double d) : v(d) {}

        friend constexpr Narrow operator+(Narrow a, Narrow b) { return Narrow(a.v + b.v); }
        friend constexpr Narrow operator-(Narrow a, Narrow b) { return Narrow(a.v - b.v); }
        friend constexpr Narrow operator*(Narrow a, Narrow b) { return Narrow(a.v * b.v); }
        friend constexpr Narrow operator/(Narrow a, Narrow b) { return Narrow(a.v / b.v); }
        friend constexpr Narrow operator-(Narrow a) { return Narrow(-a.v); }
        friend constexpr auto operator<=>(Narrow, Narrow) = default;
    };

    template<typename T>
    constexpr T narrow_square(T const &x) { return x * x; }

    template<typename M, int N, typename T>
    constexpr bool accepted = requires { derivative<M, N>(narrow_square<T>, T(1)); };

    template<typename M, int N, typename T>
    constexpr bool accepted_with_step = requires { derivative<M, N>(narrow_square<T>, T(1), T(1)); };

} // namespace

template<int D>
struct std::numeric_limits<Narrow<D>>
{
    static constexpr bool is_specialized = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_exact = false;
    static constexpr int radix = 2;
    static constexpr int digits = D;
    static constexpr Narrow<D> epsilon() { return Narrow<D>(1e-3); }
    static constexpr Narrow<D> max() { return Narrow<D>(65504.0); }
    static constexpr Narrow<D> denorm_min() { return Narrow<D>(6e-8); }
};

// A rule whose integer weights T cannot hold exactly is rejected by the constraint detail::exact_weights, so that a
// requires-expression sees it (U-10). Largest integers: Central<2> 2 (f' to f''') and 6 (f''''), Central<4> up to 56,
// Central<10> 2520 for f', Forward<6> 450 for f' and 9490 for f''; Ridders uses Central<2>'s table.
static_assert(accepted<Central<2>, 1, Narrow<11>> && accepted<Central<2>, 4, Narrow<11>>);
static_assert(accepted<Central<4>, 4, Narrow<11>> && accepted<nxx::deriv::Ridders, 4, Narrow<11>>);
static_assert(!accepted<Central<10>, 1, Narrow<11>> && !accepted_with_step<Central<10>, 1, Narrow<11>>);
static_assert(accepted<Forward<6>, 1, Narrow<11>> && !accepted<Forward<6>, 2, Narrow<11>>);
static_assert(!accepted<Backward<6>, 2, Narrow<11>>);
static_assert(accepted<nxx::deriv::Ridders, 1, Narrow<2>> && accepted<nxx::deriv::Ridders, 3, Narrow<2>>);
static_assert(!accepted<nxx::deriv::Ridders, 4, Narrow<2>> && !accepted<Central<2>, 4, Narrow<2>>);
static_assert(accepted<Central<2>, 3, Narrow<2>>);
// The one-sided extrapolations use Forward<1>'s table: 1, 2, 3 and 6 for f' to f''''.
static_assert(accepted<nxx::deriv::ForwardRidders, 3, Narrow<2>>);
static_assert(!accepted<nxx::deriv::ForwardRidders, 4, Narrow<2>>);
static_assert(accepted<nxx::deriv::BackwardRidders, 4, Narrow<11>>);
// The built-in types hold every table: Central<10> f'''' needs 6,222,216, below float's 2^24.
static_assert(accepted<Central<10>, 4, double> && accepted<Forward<6>, 4, double>);
static_assert(accepted_with_step<Central<10>, 3, float> && accepted_with_step<Backward<6>, 3, float>);
// A narrow type differentiates where its digits suffice: x^2 at 3 with h = 1/2 is exact.
static_assert(derivative(narrow_square<Narrow<11>>, Narrow<11>(3), Narrow<11>(0.5))->value.v == 6.0);

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

// The one-sided extrapolations on the user types, and BackwardRidders through its mirrored callable, which passes fn
// a T const lvalue too. x^2 at 3 with h = 1/2: the forward rule's error is linear in h for a quadratic, so the first
// column is exact, 6; f'' is exact at every level, 2.
static_assert(derivative<nxx::deriv::ForwardRidders>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
static_assert(derivative<nxx::deriv::BackwardRidders>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
static_assert(derivative<nxx::deriv::BackwardRidders, 2>(&Real::square, Real(3), Real(0.5))->value.v == 2.0);
static_assert(
    derivative<nxx::deriv::ForwardRidders>([](NoDefault const &x) { return x * x; }, NoDefault(3), NoDefault(0.5))
        ->value.v
    == 6.0);
static_assert(derivative<nxx::deriv::BackwardRidders>(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
static_assert(derivative<nxx::deriv::BackwardRidders>(LvalueOnly{}, 1.0, 0.5) == 2.0);

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

    // The one-sided extrapolations above.
    CHECK(derivative<nxx::deriv::ForwardRidders>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
    CHECK(derivative<nxx::deriv::BackwardRidders>(&Real::square, Real(3), Real(0.5))->value.v == 6.0);
    CHECK(derivative<nxx::deriv::BackwardRidders, 2>(&Real::square, Real(3), Real(0.5))->value.v == 2.0);
    CHECK(derivative<nxx::deriv::BackwardRidders>(LvalueOnly{}, 1.0, 0.5) == 2.0);
    CHECK(derivative<nxx::deriv::BackwardRidders>(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);

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

TEST_CASE("deriv: a user type with numeric_limits written member by member")
{
    // The compile-time cases above, at run time, and the default step, whose bound is the one argued for Real.
    CHECK(derivative(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
    CHECK(derivative<Forward<2>>(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
    CHECK(derivative<nxx::deriv::Ridders>(plain_square, Plain(3), Plain(0.5))->value.v == 6.0);
    CHECK(derivative<Central<2>, 2>(plain_square, Plain(3), Plain(0.5))->value.v == 2.0);
    CHECK(derivative<Central<4>, 3>(plain_square, Plain(3), Plain(0.5))->value.v == 0.0);
    auto const d = derivative(plain_square, Plain(3));
    REQUIRE(d.has_value());
    CHECK(d->value.v == doctest::Approx(6.0).epsilon(2e-10));
}
