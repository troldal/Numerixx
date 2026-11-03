// TanhSinh on infinite intervals in a type whose max() is small against 1 / eps (Phase 4's independent review,
// docs/quadrature/PHASE4_API.md), apart from test_quadrature_infinite.cpp because it compiles slowly: TanhSinh builds
// its levels for a class type without a compile-time table in constant evaluation, about a minute of compilation for
// this type, which in its own file compiles beside the rest instead of after it.

#include "quadrature_check.hpp"

#include <doctest/doctest.h>

#include <cstdlib>
#include <expected>
#include <limits>
#include <optional>

using namespace quadrature_check;

// The independent review's cases (Codex, PHASE4_API.md): a type whose max() is small against 1 / eps, here 2^40 with
// 53 bits, so that a + 1 overflows on [max(), inf) and 1 / eps does not fit. As a Boost type at run time, and as a
// double wrapper in constant evaluation, whose overflow beyond 2^40 fails the constant expression.
namespace {
    using capped_boost = mp::number<mp::cpp_bin_float<53, mp::digit_base_2, void, int, -1022, 39>>;

    // Reached only on an overflow, which the library must never cause: not constexpr, so constant evaluation fails.
    [[noreturn]] void capped_overflow()
    { std::abort(); }

    struct Capped
    {
        double v;
        static constexpr double top = 1099511627776.0; // 2^40
        constexpr Capped(int i) : v(static_cast<double>(i)) {}
        constexpr Capped(double x) : v(x) {}
        static constexpr Capped checked(double x)
        {
            if (x > top || x < -top) capped_overflow();
            return Capped(x);
        }
        friend constexpr Capped operator+(Capped x, Capped y) { return checked(x.v + y.v); }
        friend constexpr Capped operator-(Capped x, Capped y) { return checked(x.v - y.v); }
        friend constexpr Capped operator*(Capped x, Capped y) { return checked(x.v * y.v); }
        friend constexpr Capped operator/(Capped x, Capped y) { return checked(x.v / y.v); }
        friend constexpr Capped operator-(Capped x) { return Capped(-x.v); }
        friend constexpr bool operator==(Capped, Capped) = default;
        friend constexpr auto operator<=>(Capped x, Capped y) { return x.v <=> y.v; }
    };

} // namespace

template<>
struct std::numeric_limits<Capped> : std::numeric_limits<double>
{
    static constexpr Capped max() noexcept { return Capped(Capped::top); }
    static constexpr Capped infinity() noexcept { return Capped(std::numeric_limits<double>::infinity()); }
};

namespace {
    constexpr Capped capped_inf = std::numeric_limits<Capped>::infinity();
    constexpr auto capped_zero = [](Capped const &) { return Capped(0); };
    // 1/(1 + x)^2 as (1/(1 + x))/(1 + x): (1 + x)^2 itself would overflow the type near its max().
    constexpr auto capped_square = [](Capped const &x) { return Capped(Capped(1) / (Capped(1) + x)) / (Capped(1) + x); };
    // 1/eps did not fit: fn = 0 on [0, inf) failed to compile; now 0 with its floor.
    constexpr auto capped_zero_result = qd::integrate<qd::TanhSinh>(capped_zero, Capped(0), capped_inf);
    constexpr auto capped_square_result = qd::integrate<qd::TanhSinh>(capped_square, Capped(0), capped_inf);
    // a + 1 overflowed on [max(), inf): now no point inside, not_converged without a best, fn not called.
    constexpr auto capped_top_result = qd::integrate<qd::TanhSinh>(capped_square, std::numeric_limits<Capped>::max(),
                                                                   capped_inf);
} // namespace
static_assert(nxx::real<Capped> && nxx::real<capped_boost>);
static_assert(capped_zero_result && capped_zero_result->value.v == 0.0);
static_assert(capped_square_result
              && nxx::math::abs(capped_square_result->value.v - 1.0) <= capped_square_result->abs_error.v);
static_assert(!capped_top_result && capped_top_result.error().reason == qd::error::not_converged
              && !capped_top_result.error().best);

TEST_CASE("quadrature: infinite: a type whose max() is small against 1 / eps")
{
    // The run-time twins of the constant expressions above.
    CHECK(qd::integrate<qd::TanhSinh>(capped_zero, Capped(0), capped_inf) == capped_zero_result);
    CHECK(qd::integrate<qd::TanhSinh>(capped_square, Capped(0), capped_inf) == capped_square_result);
    CHECK(qd::integrate<qd::TanhSinh>(capped_square, std::numeric_limits<Capped>::max(), capped_inf)
          == capped_top_result);
    // The Boost type: [max(), inf) and (-inf, -max()] sample nothing (fn was called at inf before), and an
    // ordinary integral, 1/(1 + x)^2 on [0, inf) = 1, is honest.
    capped_boost const top = (std::numeric_limits<capped_boost>::max)();
    capped_boost const inf = std::numeric_limits<capped_boost>::infinity();
    int calls = 0;
    auto const counted = [&calls](capped_boost const &x) {
        ++calls;
        return capped_boost(capped_boost(1 / (1 + x)) / (1 + x));
    };
    qd::Failure<capped_boost> const none{ qd::error::not_converged, std::nullopt };
    CHECK(qd::integrate<qd::TanhSinh>(counted, top, inf) == std::unexpected(none));
    CHECK(qd::integrate<qd::TanhSinh>(counted, -inf, capped_boost(-top)) == std::unexpected(none));
    CHECK(calls == 0);
    auto const one = qd::integrate<qd::TanhSinh>(counted, capped_boost(0), inf);
    REQUIRE(one.has_value());
    CHECK(error_of(one->value, "1") <= one->abs_error);
}
