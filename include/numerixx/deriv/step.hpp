// Step specifications (DESIGN §7.1, D32): scalar-agnostic values, resolved per stencil and per x inside diff.
//   optimal{}          factor = eps^(1/(Order + Accuracy)) as an exact power of two
//   relative{f[, t]}   h = f * max(|x|, t); without t, h = f * |x|, and f at x == 0
//   absolute{h}        h itself
// The step is then made exactly representable: h = (x + h) - x.
//
// Deferred to phase 2: noise{eps_f} for functions computed with limited precision.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/refined.hpp>
#include <numerixx/core/scalar.hpp>

#include <algorithm>
#include <optional>

NXX_BEGIN_HEADER

namespace nxx::deriv
{
    struct optimal
    {
    };

    template<real P>
    class relative
    {
        // No std::optional member: MSVC 19.51 cannot evaluate a consteval constructor that initialises one when the
        // object is a function argument (relative{1e-5, 1.0} as a step), so the typical scale is a value plus a flag.
        tolerance<P> factor_;
        tolerance<P> typical_;
        bool         has_typical_;

    public:
        constexpr relative(tolerance<P> factor, std::optional<tolerance<P>> typical = std::nullopt) noexcept
            : factor_(factor),
              typical_(typical.value_or(factor)),
              has_typical_(typical.has_value())
        {}

        consteval relative(P factor, P typical)
            : factor_(nxx::detail::trust_me {}, factor),
              typical_(nxx::detail::trust_me {}, typical),
              has_typical_(true)
        {
            if (!tag::positive_tolerance::check(factor) || !tag::positive_tolerance::check(typical))
                nxx::detail::literal_violates_invariant("relative{factor, typical} needs finite factor > 0 and typical > 0");
        }

        constexpr tolerance<P>                factor() const noexcept { return factor_; }
        constexpr std::optional<tolerance<P>> typical() const noexcept
        { return has_typical_ ? std::optional<tolerance<P>>(typical_) : std::nullopt; }
    };

    template<real P>
    relative(P) -> relative<P>;
    template<real P>
    relative(P, P) -> relative<P>;

    template<real P>
    class absolute
    {
        tolerance<P> h_;

    public:
        constexpr absolute(tolerance<P> h) noexcept : h_(h) {}
        constexpr tolerance<P> h() const noexcept { return h_; }
    };

    template<real P>
    absolute(P) -> absolute<P>;

    namespace detail
    {
        template<real T>
        constexpr T scale_of(const T& x) noexcept
        { return x == T(0) ? T(1) : math::abs(x); }

        template<int O, int A, real T>
        constexpr T raw_step(const optimal&, const T& x) noexcept
        { return math::root_eps<T>(1, O + A) * detail::scale_of(x); }

        template<int O, int A, real T, class P>
        constexpr T raw_step(const relative<P>& r, const T& x) noexcept
        {
            const T scale = r.typical() ? (std::max)(math::abs(x), T(r.typical()->value())) : detail::scale_of(x);
            return T(r.factor().value()) * scale;
        }

        template<int O, int A, real T, class P>
        constexpr T raw_step(const absolute<P>& a, const T&) noexcept
        { return T(a.h().value()); }

        // The step actually taken: exactly representable, so x + h - x == h.
        template<int O, int A, class H, real T>
        constexpr T resolve(const H& h, const T& x) noexcept
        {
            const T raw = detail::raw_step<O, A>(h, x);
            const T xh  = x + raw;
            return xh - x;
        }
    }    // namespace detail
}    // namespace nxx::deriv

NXX_END_HEADER
