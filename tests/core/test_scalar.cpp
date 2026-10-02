// Scalars (DESIGN §6.1): a real type must have a specialised std::numeric_limits, because the default criteria, the
// power-of-two constants and math::isfinite read digits, epsilon, max and lowest from it. A type that specialises only
// scalar_traits is rejected at compile time rather than run with zero-valued limits (on x^2 - 2 the default criteria
// then stopped after at most one iteration, at errors of 0.41 to 0.59, and still reported stop_reason::criterion).
#include <numerixx/core.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <compare>
#include <limits>

namespace
{
    // A minimal real wrapper over double. Its arithmetic is enough for nxx::real; only its limits differ below.
    template<int Tag>
    struct wrapped
    {
        double v            = 0.0;
        constexpr wrapped() = default;
        constexpr wrapped(double x) noexcept : v(x) {}
        constexpr wrapped(int x) noexcept : v(x) {}
        friend constexpr wrapped operator+(wrapped a, wrapped b) noexcept { return a.v + b.v; }
        friend constexpr wrapped operator-(wrapped a, wrapped b) noexcept { return a.v - b.v; }
        friend constexpr wrapped operator*(wrapped a, wrapped b) noexcept { return a.v * b.v; }
        friend constexpr wrapped operator/(wrapped a, wrapped b) noexcept { return a.v / b.v; }
        friend constexpr wrapped operator-(wrapped a) noexcept { return -a.v; }
        constexpr wrapped&       operator*=(wrapped b) noexcept
        {
            v *= b.v;
            return *this;
        }
        constexpr wrapped& operator/=(wrapped b) noexcept
        {
            v /= b.v;
            return *this;
        }
        friend constexpr bool operator==(wrapped, wrapped) = default;
        friend constexpr auto operator<=>(wrapped a, wrapped b) noexcept { return a.v <=> b.v; }
        friend wrapped        abs(wrapped a) noexcept { return std::abs(a.v); }    // found by ADL from nxx::math
        friend bool           isfinite(wrapped a) noexcept { return std::isfinite(a.v); }
    };

    using traits_only = wrapped<0>;    // scalar_traits specialised, std::numeric_limits not
    using with_limits = wrapped<1>;    // std::numeric_limits specialised
}    // namespace

template<>
struct nxx::scalar_traits<traits_only>
{
    static constexpr bool        is_real = true;
    static constexpr traits_only epsilon() noexcept { return traits_only(0x1p-52); }
};

template<>
struct std::numeric_limits<with_limits>
{
    static constexpr bool        is_specialized = true;
    static constexpr bool        is_integer     = false;
    static constexpr bool        is_exact       = false;
    static constexpr int         digits         = std::numeric_limits<double>::digits;
    static constexpr with_limits epsilon() noexcept { return std::numeric_limits<double>::epsilon(); }
    static constexpr with_limits min() noexcept { return (std::numeric_limits<double>::min)(); }
    static constexpr with_limits max() noexcept { return (std::numeric_limits<double>::max)(); }
    static constexpr with_limits lowest() noexcept { return std::numeric_limits<double>::lowest(); }
    static constexpr with_limits infinity() noexcept { return std::numeric_limits<double>::infinity(); }
};

TEST_SUITE("core")
{
    TEST_CASE("scalars: a type without std::numeric_limits is not real, even with a scalar_traits specialisation")
    {
        static_assert(!std::numeric_limits<traits_only>::is_specialized);
        static_assert(!nxx::is_real_v<traits_only>);
        static_assert(!nxx::real<traits_only>);
        static_assert(!nxx::is_real_v<const traits_only&>);
    }

    TEST_CASE("scalars: arrays and functions are not real, and asking is not a hard error")
    {
        // numeric_limits<double[2]> is ill-formed; the trait must say false before it is instantiated, or a deleted
        // overload that should report "two ends of a real type" for a 2-D array would fail inside <limits> instead.
        static_assert(!nxx::is_real_v<double[2]>);
        static_assert(!nxx::is_real_v<const double (&)[2]>);
        static_assert(!nxx::is_real_v<double()>);
        static_assert(!nxx::real<double[2]>);
    }

    TEST_CASE("scalars: a type with a specialised std::numeric_limits is real and gets the defaults of its precision")
    {
        static_assert(nxx::is_real_v<with_limits>);
        static_assert(nxx::real<with_limits>);
        static_assert(nxx::math::root_eps<with_limits>(1, 2) == with_limits(0x1p-26));
        static_assert(nxx::floored_width {}.factor<with_limits>() == with_limits(0x1p-50));
        static_assert(nxx::math::isfinite(with_limits(3.0)));
        static_assert(nxx::bracket { with_limits(0.0), with_limits(3.0) }.hi() == with_limits(3.0));
        CHECK(nxx::bracket<with_limits>::make(with_limits(3.0), with_limits(0.0)) == nxx::bracket { with_limits(0.0), with_limits(3.0) });
    }
}
